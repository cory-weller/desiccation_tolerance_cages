#!/usr/bin/env bash

wget -L https://virginia.box.com/shared/static/hne9mcpr3nd3yur11aeg9zwfc7m5cytz.tar -O low_coverage_dna.tar && \
tar -xf low_coverage_dna.tar -d ./_dna/ && \
rm low_coverage_dna.tar
