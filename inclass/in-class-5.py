# In class 5 Question 1
def hashtag(aa):
    """
    This function takes a nucleotide trimer and returns it numeric value 0-63
    e.g. AAA return 0, AAC returns 1
    :param aa: a nucleotide trimer string
    :return: a numeric value of the trimer between 0 and 63
    """
    # Defining an empty list
    int_vec = [None]*3
    # Defining the place in the empty list
    place = 0
    # for loop that converts each element of the trimer to an integer
    for nt in aa:
        if nt == "A":
            int_vec[place] = 0
        elif nt == "C":
            int_vec[place] = 1
        elif nt == "G":
            int_vec[place] = 2
        else:
            int_vec[place] = 3

        place = place+1

    # Calculate the integer value of the trimer
    output = int_vec[0]*(4**2) + int_vec[1]*(4**1) + int_vec[2]*(4**0)

    return output

hashtag("AAA") # equals 0
hashtag("AAC") # equals 1
hashtag("ATC") # equals 13 - from example
hashtag("TTT") # last trimer - equals 63


# Question 2
# create all possible nucleotides
import itertools as it

#defining alphabet of nucleotides
alphabet = ["A", "C", "G", "T"]
# creating a vector of all possible codons
aa = [''.join(i) for i in it.product(alphabet, repeat = 3)]

def trimer_count(seq, aa):
    """

    :param seq:
    :param aa:
    :return:
    """
    loop_len = len(seq)-2
    count_vec = [0]*64
    for i in range(loop_len):
        trimer = seq[i:i+3]
        idx = hashtag(trimer)
        count_vec[idx] = count_vec[idx] + 1


    for i in range(len(count_vec)):
        if count_vec[i] > 0:
            print aa[i] + ": " + str(count_vec[i])

    return None

test_seq = "ATTATTGC"
trimer_count(test_seq, aa)
# ATT: 2
# TAT: 1
# TGC: 1
# TTA: 1
# TTG: 1