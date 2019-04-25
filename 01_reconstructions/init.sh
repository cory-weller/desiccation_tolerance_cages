#!/usr/bin/env bash

# Retrieve and unzip harp, htslib, GATK, etc
wget -L https://virginia.box.com/shared/static/ntqgqmkl5xtw10a3awvyoazptkbhyfnx.gz -O etc.tar.gz
tar -zxvf etc.tar.gz && rm etc.tar.gz
