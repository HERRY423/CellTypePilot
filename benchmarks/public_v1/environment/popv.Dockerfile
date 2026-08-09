# syntax=docker/dockerfile:1.7
FROM python:3.12.13-slim@sha256:229a2c5bfa27522db7815ea81f9bed70af17ccb9de9fc7ad142b1877b5830d36

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=1000 \
    PIP_RETRIES=10 \
    DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ \
    && rm -rf /var/lib/apt/lists/*

RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install \
      --index-url https://download.pytorch.org/whl/cpu \
      torch==2.13.0+cpu

RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install popv==0.6.1 \
    && python -m pip check \
    && python -c "import popv; print(popv.__version__)"

LABEL org.opencontainers.image.title="CellTypePilot public benchmark popV runtime" \
      org.opencontainers.image.version="0.6.1" \
      org.opencontainers.image.source="benchmarks/public_v1/environment/popv.Dockerfile"
