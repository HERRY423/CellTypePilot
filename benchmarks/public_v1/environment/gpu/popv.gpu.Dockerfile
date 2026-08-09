# syntax=docker/dockerfile:1.7
# GPU popV retrain runtime for CellTypePilot public benchmark (frozen GPU batch).
# Resolve and pin base digest on the build host, then record built_image_digest_pin
# in stack_pins.json / image_identity.json. Do not mix outputs with CPU popV.

ARG CUDA_BASE=nvidia/cuda:12.4.1-cudnn-runtime-ubuntu22.04
FROM ${CUDA_BASE}

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=1000 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-venv \
        g++ \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/bin/python3 /usr/local/bin/python \
    && ln -sf /usr/bin/pip3 /usr/local/bin/pip

# Pin torch CUDA wheel channel (cu124). Bump only with a new batch id.
RUN python -m pip install --upgrade pip \
    && python -m pip install \
        --index-url https://download.pytorch.org/whl/cu124 \
        torch==2.4.1+cu124

RUN python -m pip install popv==0.6.1 \
    && python -m pip check \
    && python -c "import torch; import popv; print('torch', torch.__version__, 'cuda', torch.cuda.is_available()); print('popv', popv.__version__)"

LABEL org.opencontainers.image.title="CellTypePilot public benchmark popV GPU runtime" \
      org.opencontainers.image.version="0.6.1-cu124" \
      org.opencontainers.image.source="benchmarks/public_v1/environment/gpu/popv.gpu.Dockerfile" \
      celltypepilot.batch="gpu_popv_retrain_v1" \
      celltypepilot.device="cuda"
