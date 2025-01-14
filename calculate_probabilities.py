import pickle
from mtokenize import tokenize

def read_tokens(fname):
    tokens = []
    with open(fname) as inp:
        for line in inp:
            tokens.append(line.rstrip().split()[0])
    return tokens

SOS, EOS = range(2)

if __name__ == "__main__":
    # The probabilities are conditioned on a prefix of length PREFIX_LENGTH
    # Why 8? No particular reason, but it's generally long enough to see the other side of a ring system.
    PREFIX_LENGTH = 8

    allowed_tokens = read_tokens("tokens.curated.txt")
    lookup_tokens = dict((x, i) for (i, x) in enumerate(allowed_tokens))

    N = len(allowed_tokens)
    dims = [N] * (PREFIX_LENGTH + 1)
    counts = {}

    # The input files were created by taking all of ChEMBL parent SMILES and
    # generating three files containing randomly ordered SMILES.
    for epoch in range(3):
        with open(f"../nbu/parent_canonical_smiles.randomize.{epoch}.smi") as inp:
            lineno = 0
            for line in inp:
                lineno += 1
                smi = line.rstrip()
                vals = [SOS] * PREFIX_LENGTH
                try:
                    for token in tokenize(smi):
                        vals.append(lookup_tokens[token])
                except KeyError:
                    continue # only allowed tokens
                vals += [EOS]

                for i, val in enumerate(vals[:-PREFIX_LENGTH]):
                    idx = vals[i+8] + (val<<6) + (vals[i+1]<<12) + (vals[i+2]<<18) + (vals[i+3]<<24) + (vals[i+4]<<30) + (vals[i+5]<<36) + (vals[i+6]<<42) + (vals[i+7]<<48)
                    counts[idx] = counts.get(idx, 0) + 1

    with open(f"../nbu/probabilities.{PREFIX_LENGTH}.pkl", "wb") as out:
        pickle.dump(obj=counts, file=out)
