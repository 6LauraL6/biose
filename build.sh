#!/bin/bash

if ! command -v docker &> /dev/null; then
    sudo apt update
    sudo apt install -y docker.io
fi

docker build -t mamorosdev/m14-uf2-bioseq-2 .