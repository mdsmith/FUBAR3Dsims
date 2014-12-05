#! /usr/bin/env python3

import numpy as NP
import argparse
import math
import treegen

def make_gen(mean, stdev):
  def omegen():
    sample = -1
    while sample <= 0:
      sample = NP.random.normal(mean, stdev)
    return sample
  return omegen

def list_to_hylist(name, pylist):
  strpylist = [str(a) for a in pylist]
  return name + " = {{" + ', '.join(strpylist) + "}};\n"

def generate_settings(num_taxa,
                      distribution,
                      out_name,
                      nonsyn,
                      syn):
    #print(num_taxa)
    #print(distribution)
    #print(out_name)
    this_set = {}
    out = open(out_name, 'w')
    buffer = []

    #treeDef = "Tree bsrel_tree = ((1,2)Interior,3,4);\n"
    treeDef = "Tree bsrel_tree = "  + treegen.get_unrooted_tree_iteratively(num_taxa) + ";\n"
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
    nonsyngens = [make_gen(nonsyn[0], nonsyn[1]),
                  make_gen(nonsyn[2], nonsyn[3])]
    #syngens = [make_gen(.4, .8), make_gen(.4, .8)]
    syngens = [ make_gen(syn[0], syn[1]),
                make_gen(syn[2], syn[3])]
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
  entry += "\t\"length\" : {{" + param_list["length"][0]
  entry += "}},\n"
  entry += "\t\"omegas\" : {"
  for n,s,p in zip(param_list["nonsyns"], param_list["syns"], param_list["props"]):
    entry += "\t{ " + str(n) + ", " + str(s) + ", "  + str(p) + "}\n"
  entry += "}\n},\n"
  return entry

def generate_taxa(tax_name):
  entry = {}
  entry["name"] = str(tax_name)
  length = -1
  while length <= 0:
    length = NP.random.normal(0.3, .2)
  entry["length"] = [str(length)]
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

def generate_all_settings(num_taxa, dist_file, out_base, nonsyn, syn):
  reps = rep_csv(dist_file)
  for i, rep in enumerate(reps):
    generate_settings(num_taxa,
                      rep,
                      out_base + "." + str(i),
                      nonsyn,
                      syn)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('num_taxa', help="the number of taxa to simulate")
  parser.add_argument('distribution_samples',
                      help="the csv of which distribution to sample")
  parser.add_argument('out_file',
                      help="the base output file name")
  parser.add_argument('--nonsyn',
                      nargs=4,
                      default=[4, 1, 0.4, 0.2],
                      type=float,
                      help="the mean and stdevs of the two nonsynonymous \
                          distributions")
  parser.add_argument('--syn',
                      nargs=4,
                      default=[0.4, 0.1, 0.4, 0.1],
                      type=float,
                      help="the mean and stdevs of the two synonymous \
                          distributions")
  args = parser.parse_args()
  generate_all_settings(int(args.num_taxa),
                        args.distribution_samples,
                        args.out_file,
                        args.nonsyn,
                        args.syn)
