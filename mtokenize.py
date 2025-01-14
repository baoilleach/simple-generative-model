chars = set('[lr%')

def tokenize(smi):
    """
    >>> tokenize_v3("Cl11%11%(111)C[C@@H](Br)I")
    ['Cl', '1', '1', '%11', '%(111)', 'C', '[C@@H]', '(', 'Br', ')', 'I']
    """
    tokens = []
    i = 0
    N = len(smi)
    while i < N:
        x = smi[i]
        if x not in chars:
            tokens.append(x)
            i += 1
        else:
            if x == 'l':
                tokens[-1] = 'Cl'
                i += 1
            elif x == 'r':
                tokens[-1] = 'Br'
                i += 1
            elif x=='[':
                j = i+1
                while smi[j] != ']':
                    j += 1
                tokens.append(smi[i:j+1])
                i += j-i + 1
            else: # %
                if smi[i+1] == '(':
                    j = i
                    while smi[j] != ')':
                        j += 1
                    tokens.append(smi[i:j+1])
                    i += j-i + 1
                else:
                    tokens.append(smi[i:i+3])
                    i += 3
    return tokens
