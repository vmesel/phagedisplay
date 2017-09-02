b = [('1', 'atgaattactcattaaaacaattaagtgtggataaaactcctgtcatatgtgttcgtcaattatgtttaacattttacatgtacatgacaaaaaaacatacattagaacataaaggacggttagtttccacatcactgagtacattccacagaaac'), ('136', 'ctgagtacattccacagaaacttagacaaaaagaacagattggtctcaaccgtgctactgcccattatttgcttgatagtggatataatgtttttttaa')]

def list_normalizer(list, overlap_len, desired_len):
    if len(list) > 1:
        normalized_list = [x for x in list[:-2]]
        seq_complete = list[-2][1]
        seq_incomplete = list[-1][1]

        if len(seq_incomplete) >= int(overlap_len):
            if len(seq_incomplete) < int(desired_len):
                missing_len = int(desired_len) - len(seq_incomplete) + 21

                print(missing_len)
                print(seq_complete)
                print(seq_incomplete)
                seq_incomplete = seq_complete[-missing_len:-21] + seq_incomplete
                print(seq_incomplete)

                normalized_list.append((list[-2][0], list[-2][1]))
                normalized_list.append((list[-1][0], seq_incomplete))
        return normalized_list
    return list

print(list_normalizer(b, 21, 156))
