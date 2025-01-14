import pickle
import tqdm
import numpy as np
import partialsmiles as ps

np.random.seed(1)

PREVENT_FULL_INVALID_STRING = False

SOS, EOS = range(2)

TOKENS = [
  '^', '$', # i.e. SOS, EOS
  'c', 'C', '(', ')', 'O', '1', '2', '=', 'N', 'n', '3', '[C@H]',
  '[C@@H]', '4', 'F', '-', 'S', '/', 'Cl', '[nH]', 's', 'o', '5', '[C@]', '#',
  '[C@@]', '\\', '[O-]', '[N+]', 'Br', '6', 'P', '7', '8', '9']

TOKEN_IDXS = list(range(len(TOKENS)))

def find_disallowed_tokens(seq, freqs):
    disallowed = set()
    smi = "".join(TOKENS[x] for x in seq)
    for i, x in enumerate(TOKENS[2:], 2):
        if freqs[i] > 0:
            try:
                ps.ParseSmiles(smi + x, partial=True)
            except ps.Error:
                disallowed.add(i)

    if PREVENT_FULL_INVALID_STRING:
        if freqs[EOS] > 0:
            try:
                ps.ParseSmiles(smi, partial=False)
            except:
                disallowed.add(EOS)
    return disallowed

def generate(all_freqs, prefix_length, temperature, avoid_invalid):
    seq = [SOS] * prefix_length
    i = 0
    while True:
        idx = (seq[i]<<6) + (seq[i+1]<<12) + (seq[i+2]<<18) + (seq[i+3]<<24) + (seq[i+4]<<30) + (seq[i+5]<<36) + (seq[i+6]<<42) + (seq[i+7]<<48)
        freqs = [all_freqs.get(idx+j, 0) for j in TOKEN_IDXS]
        if avoid_invalid:
            disallowed_tokens = find_disallowed_tokens(seq[prefix_length:], freqs)
            freqs = [all_freqs.get(idx+j, 0) if j not in disallowed_tokens else 0 for j in TOKEN_IDXS]
            if all(freq==0 for freq in freqs):
                return ""
        probs = freqs/np.sum(freqs)
        p = np.power(probs, 1.0/temperature)
        q = p / np.sum(p)
        chosen = np.random.choice(TOKEN_IDXS, 1, p=q)
        seq.append(chosen[0])
        if chosen == 1: # EOS
            break
        i += 1
    return seq[PREFIX_LENGTH:-1]

def test(seq, all_freqs):
    i = 0
    idx = (seq[i]<<6) + (seq[i+1]<<12) + (seq[i+2]<<18) + (seq[i+3]<<24) + (seq[i+4]<<30) + (seq[i+5]<<36) + (seq[i+6]<<42) + (seq[i+7]<<48)
    freqs = [all_freqs.get(idx+j, 0) for j in TOKEN_IDXS]
    return zip(TOKENS, freqs)

if __name__ == "__main__":
    PREFIX_LENGTH = 8
    TEMPERATURE = 0.8
    AVOID_INVALID = True

    with open(f"../nbu/probabilities.{PREFIX_LENGTH}.2.pkl", "rb") as inp:
        all_freqs = pickle.load(inp)
        print("Loaded...")

    ofile = f"../nbu/out.{PREFIX_LENGTH}_{TEMPERATURE}_{AVOID_INVALID}.smi"
    with open(ofile, "w") as out:
        for i in tqdm.tqdm(range(1000)):
            tokens = generate(all_freqs, PREFIX_LENGTH, TEMPERATURE, AVOID_INVALID)
            smi = "".join(TOKENS[x] for x in tokens)
            out.write(f"{smi}\n")
