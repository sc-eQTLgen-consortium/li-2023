#!/usr/bin/env bash

vcfpath=$1
snpspath=$2
savepath=$3

ml BCFtools
bcftools view --include ID==@${snpspath} ${vcfpath} > ${savepath}