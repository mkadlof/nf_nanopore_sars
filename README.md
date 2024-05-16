NextFlow pipeline for SARS-CoV-2 Nanopore data
==============================================

This pipeline is reimplementation of [pzh_pipeline_sars_nanopore](https://github.com/michallaz/pzh_pipeline_sars_nanopore) in NextFlow.

Building
------

    DOCKER_BUILDKIT=1 docker build --target main -t nf_nanopore_sars-3.0-main:latest .

You can enter docker with:

    docker run -it --rm nf_nanopore_sars-3.0-main:latest

Testing
-------

    nf-tests test
