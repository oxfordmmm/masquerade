# Use minconda base image
FROM continuumio/miniconda3:latest

# Set the working directory within the container
WORKDIR /app

# Set up bioconda
COPY env.yml /app/env.yml
RUN conda env update -n base --file env.yml && \
    conda clean -afy

# Add Python files from repo to Docker image
COPY ./src /app/src
COPY ./pyproject.toml /app/pyproject.toml
RUN pip install .
