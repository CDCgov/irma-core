FROM ghcr.io/cdcgov/irma-core/builder:latest AS builder

SHELL ["/bin/bash", "-c"]
WORKDIR /irma-core
ARG irma_core_branch

COPY . .

RUN latest=$(git tag|tail -n1) \
    && git checkout ${irma_core_branch:-$latest} \
    && cargo build --release \
    && cargo test


FROM redhat/ubi9-micro AS base

WORKDIR /app

COPY  --from=builder /irma-core/docs /app/docs
COPY   --from=builder \
    /irma-core/target/release/irma-core \
    /irma-core/Cargo.toml \
    /irma-core/Cargo.lock \
    /irma-core/LICENSE \
    /irma-core/CHANGELOG.md \
    /irma-core/CITATION.bib \
    /irma-core/CONTRIBUTORS.md \
    /irma-core/README.md \
    /app/

ENV PATH="/app:${PATH}"
WORKDIR /data