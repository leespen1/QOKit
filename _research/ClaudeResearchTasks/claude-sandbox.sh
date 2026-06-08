#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(pwd)"
CONTAINER_HOME="/home/sandbox"

NETWORK_MODE="${NETWORK_MODE:-none}"  # set to slirp4netns when you need internet

podman run --rm -it \
  --name claude-sandbox \
  --userns=keep-id \
  --security-opt=no-new-privileges \
  --cap-drop=ALL \
  --pids-limit=512 \
  --memory=4g \
  --cpus=4 \
  --read-only \
  --network="$NETWORK_MODE" \
  --tmpfs /tmp:rw,nosuid,nodev,size=1g \
  --tmpfs "$CONTAINER_HOME":rw,nosuid,nodev,size=1g \
  -v "$REPO_DIR:/work:Z,rw" \
  -w /work \
  docker.io/library/node:20-bookworm \
  bash -lc '
    set -euo pipefail

    export HOME='"$CONTAINER_HOME"'
    export NPM_CONFIG_PREFIX="$HOME/.npm-global"
    export PATH="$NPM_CONFIG_PREFIX/bin:$PATH"

    npm config set fund false >/dev/null
    npm config set update-notifier false >/dev/null

    npm i -g @anthropic-ai/claude-code
    exec claude
  '
