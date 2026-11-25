"""Fortran 77 fixed-form -> C++ helpers for the CVRQD routines (expressions, DATA arrays)."""
import re, sys

def logical_lines(path):
    """Join continuation lines; drop comments. Returns list of (first_line_no, text)."""
    out = []
    for no, raw in enumerate(open(path), 1):
        line = raw.rstrip('\n').replace('\t', '      ')
        if not line.strip():
            continue
        if line[0] in 'cC*!':
            continue
        if line.lstrip().startswith('!'):
            continue
        # strip trailing ! comments
        if '!' in line:
            line = line[:line.index('!')]
        if len(line) > 5 and line[5] not in ' 0' and line[:5].strip() == '':
            out[-1] = (out[-1][0], out[-1][1] + line[6:].strip())
        else:
            out.append((no, line[6:].strip() if len(line) > 6 else ''))
    return out

NUM = re.compile(r'(?<![A-Za-z_0-9])(\d+\.\d*|\.\d+|\d+)([eEdD][+-]?\d+)?')

def tokens(s):
    i = 0
    toks = []
    while i < len(s):
        c = s[i]
        if c.isspace():
            i += 1; continue
        m = NUM.match(s, i)
        if m and (c.isdigit() or c == '.') :
            toks.append(('num', m.group(0))); i = m.end(); continue
        m = re.match(r'[A-Za-z_][A-Za-z_0-9]*', s[i:])
        if m:
            toks.append(('id', m.group(0).lower())); i += m.end(); continue
        if s.startswith('**', i):
            toks.append(('pow', '**')); i += 2; continue
        toks.append(('op', c)); i += 1
    return toks

def num_cpp(t):
    m = NUM.fullmatch(t)
    mant, exp = m.group(1), m.group(2)
    if '.' not in mant and not exp:
        return mant  # integer
    if exp and exp[0] in 'dD':
        return mant + 'e' + exp[1:] if True else None
    # default REAL: single precision
    return 'F(' + mant + (exp if exp else '') + ')'

def conv(s):
    """Convert a Fortran expression to C++ (x**n -> ipow(x,n), literals typed)."""
    toks = tokens(s)
    out = []  # list of strings, each a 'primary' or operator
    i = 0
    def primary_end_index():
        pass
    res = []
    # build output as list of segments; handle ** by taking the last primary
    for k, (kind, val) in enumerate(toks):
        if kind == 'num':
            res.append(num_cpp(val))
        elif kind == 'id':
            res.append(val)
        elif kind == 'pow':
            res.append('**')
        else:
            res.append(val)
    # now resolve '**': base is previous primary (id/num/ parenthesised group), exponent next primary
    def resolve(seq):
        out = []
        i = 0
        while i < len(seq):
            t = seq[i]
            if t == '**':
                # pop base
                base = out.pop()
                if base == ')':
                    depth = 1; grp = [')']
                    while depth:
                        x = out.pop(); grp.insert(0, x)
                        if x == ')': depth += 1
                        elif x == '(': depth -= 1
                    # function call name?
                    if out and re.fullmatch(r'[a-z_][a-z_0-9]*', out[-1]):
                        grp.insert(0, out.pop())
                    base = ''.join(grp)
                e = seq[i+1]
                sign = ''
                if e in '+-':
                    sign = e; i += 1; e = seq[i+1]
                assert re.fullmatch(r'\d+', e), (seq, e)
                out.append('ipow(' + base + ',' + sign + e + ')')
                i += 2
                continue
            out.append(t); i += 1
        return ''.join(out)
    r = resolve(res)
    r = r.replace('.gt.', '>').replace('.lt.', '<')
    return r

def data_values(text):
    """Values of a DATA list: 'a,b,c' possibly with r*val repeats."""
    vals = []
    for part in text.split(','):
        part = part.strip()
        if not part: continue
        if '*' in part:
            r, v = part.split('*'); vals += [v.strip()] * int(r)
        else:
            vals.append(part)
    return vals

if __name__ == '__main__':
    print(conv(sys.argv[1]))
