
def dot2ct(dot):
    """
    dot                     -- Dotbracket structure
    
    Convert dotbracket structure to list
    ..((..))..  =>  [(3, 8), (4, 7)]
    """
    
    stack = []
    ctList = []
    
    for idx,symbol in enumerate(dot):
        if symbol in ('(', '[', '{', '<'):
            stack.append( (idx+1, symbol) )
        elif symbol in '.-_=:,':
            pass
        else:
            for i in range( len(stack)-1, -1, -1 ):
                if stack[i][1]+symbol in ("()", "[]", r"{}", "<>"):
                    ctList.append( (stack[i][0], idx+1) )
                    del stack[i]
                    break
                else:
                    pass
    
    if len(stack) != 0:
        sys.stderr.writelines("Error: Bad dotbracket structure\n")
        raise NameError("Error: Bad dotbracket structure")
    
    ctList.sort()
    return ctList

def dot2bpmap(dot):
    """
    dot                 -- Dotbracket structure
    
    Convert dotbracket structure to dictionary
    ..(((....)))..   =>    {3: 12, 4: 11, 5: 10, 10: 5, 11: 4, 12: 3}
    """
    stack = []
    bpmap = {}
    
    for idx,symbol in enumerate(dot):
        if symbol in ('(', '[', '{', '<'):
            stack.append( (idx+1, symbol) )
        elif symbol in '.-_=:,': ####
            pass
        else:
            for i in range( len(stack)-1, -1, -1 ):
                if stack[i][1]+symbol in ("()", "[]", r"{}", "<>"):
                    bpmap[idx+1] = stack[i][0]
                    bpmap[stack[i][0]] = idx+1
                    del stack[i]
                    break
                else:
                    pass
    
    if len(stack) != 0:
        sys.stderr.writelines("Error: Bad dotbracket structure\n")
        raise NameError("Error: Bad dotbracket structure")
    
    return bpmap

def ct2dot(ctList, length):
    """
    ctList              -- paired-bases: [(3, 8), (4, 7)]
    length              -- Length of structure
    
    Convert ctlist structure to dot-bracket
    [(3, 8), (4, 7)]  => ..((..))..
    """
    dot = ['.']*length
    if len(ctList) == 0:
        return "".join(dot)
    ctList = sorted(ctList, key=lambda x:x[0])
    ctList = [ it for it in ctList if it[0]<it[1] ]
    pseudo_duplex = parse_pseudoknot(ctList)
    for l,r in ctList:
        dot[l-1] = '('
        dot[r-1] = ')'
    dottypes = [ '<>', r'{}', '[]' ]
    if len(pseudo_duplex)>len(dottypes):
        print("Warning: too many psudoknot type: %s>%s" % (len(pseudo_duplex),len(dottypes)))
    for i,duplex in enumerate(pseudo_duplex):
        for l,r in duplex:
            dot[l-1] = dottypes[i%3][0]
            dot[r-1] = dottypes[i%3][1]
    return "".join(dot)

def find_stem(dot, max_stem_gap=3, min_stem_len=5):
    """
    dot                 -- Dotbracket structure
    max_stem_gap        -- Maximun interior loop size
    min_stem_len        -- Minimun stem length
    
    Find stem loop structure
    """
    ctList = dot2ct(dot)
    
    vs = []
    ls = 0; le = 0; rs = 0; re = 0
    
    for pair in ctList:
        if ls == 0:
            ls = le = pair[0]
            rs = re = pair[1]
        elif abs(pair[0]-le)-1>max_stem_gap or abs(rs-pair[1])-1>max_stem_gap:
            if le-ls+1>=min_stem_len and re-rs+1>=min_stem_len:
                vs.append( (ls, le, rs, re) )
            ls = le = pair[0]
            rs = re = pair[1]
        else:
            le = pair[0]
            rs = pair[1]
    if le-ls+1>=min_stem_len and re-rs+1>=min_stem_len:
        vs.append( (ls, le, rs, re) )
    return vs


def trim_stem(dot, stem, min_fix_stem_len=2):
    """
    dot                 -- Dotbracket structure
    stem                -- [left_start, left_end, right_start, right_end]
    min_fix_stem_len    -- Minimun stem length in the stem loop boundary
    
    Trim short stem from stem loop boundary
    """
    subDot = dot[ stem[0]-1:stem[3] ]
    fix_stems = find_stem(subDot, max_stem_gap=0, min_stem_len=1)
    #print fix_stems
    
    fix_stems.sort(key=lambda x: x[0])
    for i in range(len(fix_stems)):
        cstem = fix_stems[i]
        if cstem[1]-cstem[0]+1>=min_fix_stem_len:
            return (cstem[0]+stem[0]-1, stem[1], stem[2], cstem[3]+stem[0]-1)
    
    return None

def find_stem_loop(dot, max_loop_len=4, max_stem_gap=3, min_stem_len=5):
    """
    dot                 -- Dotbracket structure
    max_loop_len        -- Maximun loop size
    max_stem_gap        -- Maximun interior loop size
    min_stem_len        -- Minimun stem length
    
    Find stem loops structure
    """
    stems = find_stem(dot, max_stem_gap=max_stem_gap, min_stem_len=min_stem_len)
    
    stem_loops = []
    for stem in stems:
        loopLen = stem[2]-stem[1]-1
        if loopLen <= max_loop_len:
            stem_loops.append( stem )
    
    return stem_loops


def find_bulge_interiorLoop(dot, stem):
    """
    dot                 -- Dotbracket structure
    stem                -- [left_start, left_end, right_start, right_end]
    
    Fint bulge and interiorLoop from stem loop

    Return bulge_list, interiorLoop_list
    
    Example:
        seq= "UCUAGGUGAUUUCUGUGAAAUCGAGCCCACUUGAUUGUUUCUGUGAAACACUCUA"
        dot = "....(((((((((...))))))..))).....((.((((((...)))))).)).."
        find_stem(dot, max_stem_gap=1, min_stem_len=5)
        stemLoop = find_stem_loop(dot, max_loop_len=3, max_stem_gap=3, min_stem_len=5)
        bulge, interiorLoop = find_bulge_interiorLoop(dot, stemLoop[0])
    """
    bulge = []
    interiorLoop = []
    
    ls,le,rs,re = stem
    i = le
    j = rs
    while i >= ls and j <= re:
        while i>=ls and j<=re and dot[i-1]!='.' and dot[j-1]!='.':
            i -= 1
            j += 1
        if i < ls or j > re:
            break
        
        if dot[i-1]!='.' and dot[j-1]=='.':
            start = j
            j += 1
            while j<=re and dot[j-1]=='.':
                j += 1
            bulge.append( (start, j-1) )
        
        if dot[i-1]=='.' and dot[j-1]!='.':
            start = i
            i -= 1
            while i>=ls and dot[i-1]=='.':
                i -= 1
            bulge.append( (i+1, start) )
        
        if dot[i-1]=='.' and dot[j-1]=='.':
            l_s = i
            r_s = j
            i -= 1
            j += 1
            while j<=re and dot[j-1]=='.':
                j += 1
            while i>=ls and dot[i-1]=='.':
                i -= 1
            interiorLoop.append( (i+1,l_s,r_s,j-1) )
    
    return bulge, interiorLoop
