import sys
argv = sys.argv

tag_ids = argv[1].split(',')
tag_files = argv[2].split(',')
id_out = argv[3]

id_len = len(tag_ids)
file_len = len(tag_files)

if id_len != file_len:
    raise ValueError("Arrays must have the same size")

i=0
with open('%s_samples.txt' % id_out,'w') as f:
    while i < id_len:
        id = tag_ids[i]
        file = tag_files[i]
        f.write('%s\t%s\n' % (id, file))
        i += 1
