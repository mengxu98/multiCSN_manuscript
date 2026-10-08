#!/usr/bin/env python3

import importlib
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from itertools import permutations


def run_inference(expression, out_file, assoc="correlation", pos_density=0.01, neg_density=0.01):
    repo_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(repo_dir))

    if "rpy2" not in sys.modules:
        rpy2 = types.ModuleType("rpy2")
        robjects = types.ModuleType("rpy2.robjects")
        packages = types.ModuleType("rpy2.robjects.packages")

        class DummyPackage:
            def __getattr__(self, name):
                raise RuntimeError(f"Optional rpy2-backed functionality is unavailable: {name}")

        def importr(_name):
            return DummyPackage()

        packages.importr = importr
        robjects.packages = packages
        robjects.numpy2ri = types.SimpleNamespace(activate=lambda: None)
        rpy2.robjects = robjects
        sys.modules["rpy2"] = rpy2
        sys.modules["rpy2.robjects"] = robjects
        sys.modules["rpy2.robjects.packages"] = packages

    correlation = importlib.import_module("pysrc.associations.correlation")
    dotprod = importlib.import_module("pysrc.associations.dotprod")
    proprho = importlib.import_module("pysrc.associations.proprho")
    signed = importlib.import_module("pysrc.graphlearning.signed")
    unsigned = importlib.import_module("pysrc.graphlearning.unsigned")

    def find_bs_upper_bound(k, d, density):
      value = 5
      while True:
        w = unsigned.learn_ladmm(k, d, value, value)
        density_est = np.count_nonzero(w) / len(w)
        if density_est > density:
          return value
        value *= 1.5

    def binary_search(k, d, density_pos, density_neg):
      apos_done = density_pos <= 0
      aneg_done = density_neg <= 0
      apos_min = 0
      aneg_min = 0
      apos_max = find_bs_upper_bound(k, d, density_pos) if density_pos > 0 else 0
      aneg_max = find_bs_upper_bound(-k, d, density_neg) if density_neg > 0 else 0

      wpos = np.zeros((len(k), 1))
      wneg = np.zeros((len(k), 1))
      for _ in range(50):
        apos = (apos_min + apos_max) / 2
        aneg = (aneg_min + aneg_max) / 2
        if aneg == 0:
          wpos = unsigned.learn(k, d, apos)
        elif apos == 0:
          wneg = unsigned.learn(-k, d, aneg)
        else:
          wpos, wneg = signed.learn(k, d, apos, aneg, lpos_init="zeros", lneg_init="zeros")

        pos_est = np.count_nonzero(wpos) / len(wpos)
        neg_est = np.count_nonzero(wneg) / len(wneg)

        if not apos_done:
          if abs(density_pos - pos_est) < 1e-2:
            apos_done = True
          elif density_pos < pos_est:
            apos_max = apos
          else:
            apos_min = apos
        if not aneg_done:
          if abs(density_neg - neg_est) < 1e-2:
            aneg_done = True
          elif density_neg < neg_est:
            aneg_max = aneg
          else:
            aneg_min = aneg
        if apos_done and aneg_done:
          break
      return wpos, -wneg

    def convert_df(gene_names, lpos, lneg):
      gene1 = [i for i, _ in permutations(gene_names, r=2)]
      gene2 = [j for _, j in permutations(gene_names, r=2)]
      matrix = squareform(np.squeeze(lpos + lneg))
      edge_weights = [matrix[i, j] for i, j in permutations(range(len(gene_names)), r=2) if i != j]
      result = pd.DataFrame({"Gene1": gene1, "Gene2": gene2, "EdgeWeight": edge_weights})
      return result[result["EdgeWeight"] != 0]

    assocs = {
      "dotprod": dotprod.calc,
      "correlation": correlation.calc,
      "proprho": proprho.calc,
    }

    expression_matrix = expression.to_numpy(dtype=float)
    gene_names = expression.index.to_numpy()
    assoc_fun = assocs[assoc]
    nonzero_mask = np.count_nonzero(expression_matrix, axis=1) != 0
    matrix = expression_matrix[nonzero_mask, :]
    genes = gene_names[nonzero_mask]
    kernel = assoc_fun(matrix)
    upper = kernel[np.triu_indices_from(kernel, k=1)]
    upper /= np.max(np.abs(upper))
    diag = np.diag(kernel) / np.max(np.abs(upper))
    wpos, wneg = binary_search(upper, diag, pos_density, neg_density)
    result = convert_df(genes, wpos, wneg)
    result.to_csv(out_file, sep="\t", index=False)

