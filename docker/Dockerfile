FROM ubuntu:24.04

# Install necessary packages
RUN apt-get update && apt-get install -y wget curl

# Download and install pixi
RUN curl -Ls https://github.com/prefix-dev/pixi/releases/download/v0.41.4/pixi-x86_64-unknown-linux-musl -o /usr/local/bin/pixi \
    && chmod +x /usr/local/bin/pixi

# Create the /metabolite_search directory
RUN mkdir -p /metabolite_search

# Set the environment variable for Python in pixi's default environment
ENV QUARTO_PYTHON='/metabolite_search/.pixi/envs/default/bin/python'

ENTRYPOINT ["/usr/local/bin/pixi", "r"]