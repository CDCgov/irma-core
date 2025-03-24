FROM redhat/ubi8:latest

ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH

RUN yum update -y && yum install -y zip git which gcc && yum clean all

ARG rustupSha256='a3339fb004c3d0bb9862ba0bce001861fe5cbde9c10d16591eb3f39ee6cd3e7f'
ARG url='https://static.rust-lang.org/rustup/archive/1.28.1/x86_64-unknown-linux-gnu/rustup-init'

RUN curl --proto '=https' --tlsv1.2 -sSf -o rustup-init "$url" \
    && echo "${rustupSha256} *rustup-init" | sha256sum -c - \
    && chmod +x rustup-init \
    && ./rustup-init -y --no-modify-path --profile minimal --default-toolchain nightly \
    && rm rustup-init \
    && chmod -R a+w $RUSTUP_HOME $CARGO_HOME \
    && rustup --version \
    && cargo --version \
    && rustc --version
