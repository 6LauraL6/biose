#!/bin/bash
set -e

## node

if [ ! -d "~/.nvm/nvm.sh" ]; then
    curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.1/install.sh | bash
fi

source ~/.nvm/nvm.sh

if ! command -v npm &> /dev/null; then
    nvm install --lts
fi

if [ ! -d "node_modules" ]; then
     npm install
fi

## flask

if ! command -v pip &> /dev/null; then
    sudo apt update
    sudo apt install -y python3-pip python3.10-venv
fi

if [ ! -d ".venv" ]; then
    python3 -m venv .venv
    .venv/bin/pip install -r requirements-dev.txt
fi