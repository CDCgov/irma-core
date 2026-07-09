#!/bin/bash

app=irma-core

command -v "cargo" > /dev/null 2>&1 || {
    echo "ERROR: you may need to install Rust nightly!"
    exit 1
}

if [ -n "$1" ]; then
    if grep -Eq "not a recognized processor" <(rustc -C "target-cpu=$1" --print host-tuple 2>&1); then
        echo "WARNING: cannot find target-cpu '$1', skipping"
    else
        export RUSTFLAGS="${RUSTFLAGS:+$RUSTFLAGS }-C target-cpu=$1"
    fi
fi

cargo build --profile prod \
    && cp target/prod/$app . \
    || {
        echo "ERROR: Installation failed!"
        exit 1
    }
