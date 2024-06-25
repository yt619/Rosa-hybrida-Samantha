#!/usr/bin/env python

import sys
import pysam
import collections
from optparse import OptionParser

VERSION = "0.1"

def process_bed(file):
  ret = collections.OrderedDict()
  with open(file, "r") as f:
    for line in f:
      lst = line.strip().split("\t")
      name, st, en = lst[0], int(lst[1]), int(lst[2])
      if name not in ret: ret[name] = []
      ret[name].append((st, en))
  return ret

def rewrite_header(input_bam, bed_dict):
  header = input_bam.header.copy()
  header['SQ'] = []
  for key in bed_dict:
    if len(bed_dict[key]) == 1:
      tmp = {}
      tmp['SN'], tmp['LN'] = key, bed_dict[key][0][1]
      header['SQ'].append(tmp)
    else:
      for i in range(len(bed_dict[key])):
        tmp = {}
        tmp['SN'] = key + "_" + str(i + 1)
        tmp['LN'] = bed_dict[key][i][1] - bed_dict[key][i][0] + 1
        header['SQ'].append(tmp)
  return header

def reconfirm_aln(aln, order, bed):
  name = aln.reference_name
  num = len(bed[name])
  if num == 1:
    return order[name], aln.reference_start
  else:
    aln_length = aln.query_alignment_length
    for i in range(num):
      if aln.reference_start >= (bed[name][i][0] - 1) and \
      (aln.reference_start + aln_length) <= (bed[name][i][1] - 1):
        name1 = name + "_" + str(i + 1)
        return order[name1], aln.reference_start - (bed[name][i][0] - 1)
    return None, None

def resplit_aln(aln1, aln2, order, bed):
  (id1, st1) = reconfirm_aln(aln1, order, bed)
  if id1 == None: return None, None
  (id2, st2) = reconfirm_aln(aln2, order, bed)
  if id2 == None: return None, None

  aln1.reference_id, aln1.reference_start = id1, st1
  aln1.next_reference_id, aln1.next_reference_start = id2, st2
  aln2.reference_id, aln2.reference_start = id2, st2
  aln2.next_reference_id, aln2.next_reference_start = id1, st1
  return aln1, aln2

def rewrite_aln(input_bam, output_bam, header, bed_dict):
  header_order = dict(map(lambda x: x[::-1], enumerate(map(lambda x: x['SN'], header['SQ']))))
  flag = 1
  for aln in input_bam:
    if flag == 1:
      aln1 = aln
      flag += 1
    else:
      aln2 = aln
      flag = 1
      (aln1, aln2) = resplit_aln(aln1, aln2, header_order, bed_dict)
      if aln1 and aln2:
        output_bam.write(aln1)
        output_bam.write(aln2)

def process_bam(options, bed_dict):
  input_bam = pysam.AlignmentFile(options.input_bam, "rb")
  header = rewrite_header(input_bam, bed_dict)
  output_bam = pysam.AlignmentFile(options.output_bam, "wb", header = header)
  rewrite_aln(input_bam, output_bam, header, bed_dict)
  input_bam.close()
  output_bam.close()

def process(options):
  if not options.bed_file:
    print("Please specify bed file!")
    return
  bed_dict = process_bed(options.bed_file)

  if not options.input_bam or not options.output_bam:
    print("Please specify input and output file!")
    return
  process_bam(options, bed_dict)

def parse_command():
  usage = "Modify bam file according to HiC data\n\npython BamModifier.py -i input.bam -o output.bam -b hic.bed"
  parser = OptionParser(usage = usage, version = VERSION)
  parser.add_option("-i", dest = "input_bam", help = "input bam file")
  parser.add_option("-o", dest = "output_bam", help = "output modified bam file")
  parser.add_option("-b", dest = "bed_file", help = "bed file including breakpoints")
  return parser.parse_args()

def main():
  (options, args) = parse_command()
  options.version = VERSION
  process(options)

if __name__ == "__main__":
  main()
