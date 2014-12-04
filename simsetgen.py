#! /usr/bin/env python3

import numpy as NP
import argparse
import math
from treegen import get_unrooted_tree

def make_gen(mean, stdev):
  def omegen():
    sample = -1
    while sample <= 0:
      sample = NP.random.normal(mean, stdev)
    return sample
  return omegen

def list_to_hylist(name, pylist):
  strpylist = [str(a) for a in pylist]
  return name + " = {" + ', '.join(strpylist) + "};\n"

def generate_settings(num_taxa,
                      distribution,
                      out_name):
    #print(num_taxa)
    #print(distribution)
    #print(out_name)
    this_set = {}
    out = open(out_name, 'w')
    buffer = []

    #treeDef = "Tree bsrel_tree = ((1,2)Interior,3,4);\n"
    treeDef = "Tree bsrel_tree = "  + get_unrooted_tree(num_taxa) + ";\n"
    codoDef = "nuc3x4 = {4,3} [\"0.25\"];\n"
    nucbDef =   "nucleotide_bias_settings = { // relative to AG\n" \
                "\"AC\": 0.25,\n" \
                "\"AT\": 0.25,\n" \
                "\"CG\": 0.4,\n" \
                "\"CT\": 2.0,\n" \
                "\"GT\": 0.3\n" \
                "};\n"
    buffer.append(treeDef)
    buffer.append(codoDef)
    buffer.append(nucbDef)

    # make the distributions from which omegas are sampled
    # closures:
    nonsyngens = [make_gen(.4, .1), make_gen(4, .25)]
    syngens = [make_gen(1, .8), make_gen(1, .8)]
    nonsynsamples = []
    synsamples = []
    for i in distribution:
      nonsynsamples.append(nonsyngens[int(i) - 1]())
      synsamples.append(syngens[int(i) - 1]())
    buffer.append(list_to_hylist("nonsynvals", nonsynsamples))
    buffer.append(list_to_hylist("synvals", synsamples))

    buffer.append("bsrel_settings = { \n")

    #for tax_name in range(1, num_taxa):
    for tax_name in range(1, (2*num_taxa) - 2):
        this_set[str(tax_name)] = generate_taxa(tax_name)

    for tax_name in range(1, (2*num_taxa) - 2):
      buffer.append(params_to_string_alpha(this_set[str(tax_name)]))

    buffer.append("};\n")

    repsDef = "replicates = 1;\n"
    seqlDef = "codons = " + str(len(distribution)) + ";\n"

    buffer.append(repsDef)
    buffer.append(seqlDef)

    out.writelines(buffer)
    out.close()
    return this_set

def params_to_string_alpha(param_list):
    entry = "\"" + param_list["name"] + "\" : {"
    entry += "\t\"omegas\" : {"
    for n,s,p in zip(param_list["nonsyns"], param_list["syns"], param_list["props"]):
        entry += "\t{ " + str(n) + ", " + str(s) + ", "  + str(p) + "}\n"
    entry += "}\n},\n"
    return entry

def generate_taxa(tax_name):
  entry = {}
  entry["name"] = str(tax_name)
  entry["nonsyns"] = [math.e]
  entry["syns"] = [math.e]
  entry["props"] = [1]
  return entry

def rep_csv(csv_file_name):
  dist_reps = []
  fh = open(csv_file_name, 'r')
  lines = fh.readlines()
  if lines:
    line1 = lines[0]
    tokens = [a.strip() for a in line1.split(',')]
    dist_reps = [[] for _ in tokens]
  for line in lines:
    tokens = [a.strip() for a in line.split(',')]
    for i,token in enumerate(tokens):
      dist_reps[i].append(token)
  return dist_reps

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('num_taxa', help="the number of taxa to simulate")
    parser.add_argument('distribution_samples',
                        help="the csv of which distribution to sample")
    parser.add_argument('out_file',
                        help="the base output file name")
    args = parser.parse_args()
    distribution_reps = rep_csv(args.distribution_samples)
    for i, distribution in enumerate(distribution_reps):
      generate_settings(int(args.num_taxa),
                        distribution,
                        args.out_file + "." + str(i))
