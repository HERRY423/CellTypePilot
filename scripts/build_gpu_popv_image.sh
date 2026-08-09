#!/usr/bin/env bash
# Build the pinned GPU popV image and record digests into environment/gpu/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GPU_ENV="${ROOT}/benchmarks/public_v1/environment/gpu"
DOCKERFILE="${GPU_ENV}/popv.gpu.Dockerfile"
IMAGE_TAG="celltypepilot-popv-gpu:0.6.1-cu124"
BASE_REF="nvidia/cuda:12.4.1-cudnn-runtime-ubuntu22.04"

if ! command -v docker >/dev/null 2>&1; then
  echo "docker not found" >&2
  exit 1
fi

echo "Pulling base image for digest pin: ${BASE_REF}"
docker pull "${BASE_REF}"
BASE_DIGEST="$(docker image inspect "${BASE_REF}" --format '{{index .RepoDigests 0}}')"
echo "Base: ${BASE_DIGEST}"

echo "Building ${IMAGE_TAG}"
docker build \
  --pull=false \
  -f "${DOCKERFILE}" \
  --build-arg "CUDA_BASE=${BASE_REF}" \
  -t "${IMAGE_TAG}" \
  "${ROOT}"

IMAGE_ID="$(docker image inspect "${IMAGE_TAG}" --format '{{.Id}}')"
REPO_DIGEST="$(docker image inspect "${IMAGE_TAG}" --format '{{index .RepoDigests 0}}' 2>/dev/null || true)"
CREATED="$(docker image inspect "${IMAGE_TAG}" --format '{{.Created}}')"

python3 - <<PY
import json
from datetime import datetime, timezone
from pathlib import Path

gpu_env = Path(r"""${GPU_ENV}""")
payload = {
    "schema_version": "celltypepilot.gpu-image-identity.v1",
    "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
    "image_tag": "${IMAGE_TAG}",
    "image_id": "${IMAGE_ID}",
    "repo_digest": "${REPO_DIGEST}" or None,
    "base_ref": "${BASE_REF}",
    "base_digest": "${BASE_DIGEST}",
    "dockerfile": "benchmarks/public_v1/environment/gpu/popv.gpu.Dockerfile",
    "created": "${CREATED}",
    "batch_id": "gpu_popv_retrain_v1",
    "note": "Recorded on build host. Push to registry then refresh repo_digest if empty.",
}
out = gpu_env / "image_identity.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(f"Wrote {out}")

pins_path = gpu_env / "stack_pins.json"
pins = json.loads(pins_path.read_text(encoding="utf-8"))
pins["container_image"]["base_image_digest_pin"] = "${BASE_DIGEST}".split("@", 1)[-1] if "@" in "${BASE_DIGEST}" else "${BASE_DIGEST}"
pins["container_image"]["built_image_digest_pin"] = "${IMAGE_ID}"
pins["container_image"]["repo_digest"] = "${REPO_DIGEST}" or None
pins["frozen_at_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
pins_path.write_text(json.dumps(pins, indent=2) + "\n", encoding="utf-8")
print(f"Updated {pins_path}")
PY

echo "Smoke GPU runtime inside container"
docker run --rm --gpus all "${IMAGE_TAG}" \
  python -c "import torch, popv; assert torch.cuda.is_available(), 'CUDA not visible'; print(torch.__version__, popv.__version__, torch.cuda.get_device_name(0))"

echo "Done. Image ${IMAGE_TAG} ready for gpu_popv_retrain_v1 workers."
