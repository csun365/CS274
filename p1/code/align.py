"""
align.py
This file implements the Needleman-Wunsch and Smith-Waterman algorithms for global and local sequence alignment, respectively. 
Author: Christopher Sun

Usage: python align.py input_file_path output_file_path
"""

import sys

#### ------ USEFUL FUNCTIONS ------- ####
# Function to determine equality between two floating point numbers
# Inputs:
#       a: a floating point number
#       b: a floating point number
# Returns: 
#       a boolean signifying equality
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    
#### ------- CLASSES ------- ####
class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    """
    # Function to initialize the MatchMatrix class
    # Inputs:
    #       None
    # Returns:
    #       None
    def __init__(self):
        self.M = {}

    # Function to set the score of an entry in the match matrix
    # Inputs:
    #       a: the first letter
    #       b: the second letter
    #       score: the score of an alignment between letters a and b
    # Returns: 
    #       None
    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """
        self.M[(a, b)] = score

    # Function to get the score of an entry in the match matrix
    # Inputs:
    #       a: the first letter
    #       b: the second letter
    # Returns:
    #       the score of an alignment between letters a and b
    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """
        return self.M[(a, b)]

class Score(object):
    """
    Object to store an entry in the Score Matrix
    """
    # Function to initialize the MatchMatrix class
    # Inputs:
    #       score: the numerical value of the score
    #       pointers: a list of pointers (useful for traceback)
    # Returns:
    #       None
    def __init__(self, score, pointers):
        self.score = score
        self.pointers = pointers
    
class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """

    # Function to initialize the ScoreMatrix class
    # Inputs:
    #       name: identifier for the score matrix - Ix, Iy, or M
    #       nrow: number of rows in the matrix
    #       ncol: number of columns in the matrix
    # Returns:
    #       None
    def __init__(self, name, nrow, ncol):
        self.name = name 
        self.nrow = nrow
        self.ncol = ncol
        self.score_matrix = [[Score(0, []) for i in range(ncol)] for j in range(nrow)]

    # Function to get the score at a particular location
    # Inputs:
    #       row: the desired row
    #       col: the desired column
    # Returns:
    #       the score at index (row, col) of the matrix
    def get_score(self, row, col):
        return self.score_matrix[row][col].score

    # Function to set the score at a particular location
    # Inputs:
    #       row: the desired row
    #       col: the desired column
    #       score: the score to be set at (row, col) of the matrix
    # Returns:
    #       None
    def set_score(self, row, col, score):    
        self.score_matrix[row][col].score = score

    # Function to get the pointers at a particular location
    # Inputs:
    #       row: the desired row
    #       col: the desired column
    # Returns:
    #       the pointers at index (row, col) of the matrix
    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to
        This should be formatted as a list of tuples:
         ex. [(1,1), (1,0)]
        """
        return self.score_matrix[row][col].pointers

    # Function to get the pointers at a particular location
    # Inputs:
    #       row: the desired row
    #       col: the desired column
    #       pointer: the pointer to be appended to the list of pointers at (row, col) of the matrix
    # Returns:
    #       None
    def set_pointers(self, row, col, pointer):
        self.score_matrix[row][col].pointers.append(pointer)

    # Function to display the matrix scores
    # Inputs:
    #       None
    # Returns:
    #       None
    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        """
        print("{}=".format(self.name))
        for i in range(self.nrow):
            for j in range(self.ncol):
                print("\t{:.2f}".format(self.score_matrix[i][j].score), end="")
            print()

    # Function to display the matrix pointers
    # Inputs:
    #       None
    # Returns:
    #       None
    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry in the score matrix. Use this for debugging!
        """
        print("Pointers for {}=".format(self.name))
        for i in range(self.nrow):
            for j in range(self.ncol):
                print("\t" + str(self.score_matrix[i][j].pointers), end="")
            print()

class AlignmentParameters(object):
    """
    Object to hold a set of alignment parameters from an input file.
    """

    # Function to define alignment parameters
    # Inputs:
    #       None
    # Returns:
    #       None
    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False 
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix()

    # Function to load user-specified alignment parameters from a file
    # Inputs:
    #       input_file: path to the input file
    # Returns:
    #       None
    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """
        with open(input_file, "r") as f:
            contents = [i.strip() for i in f.readlines()]
        self.seq_a = contents[0]
        self.seq_b = contents[1]
        self.num_rows = len(self.seq_a) + 1
        self.num_cols = len(self.seq_b) + 1
        self.global_alignment = not bool(int(contents[2]))
        self.dx, self.ex, self.dy, self.ey = tuple([float(i) for i in contents[3].split()])
        self.len_alphabet_a = int(contents[4])
        self.alphabet_a = contents[5]
        self.len_alphabet_b = int(contents[6])
        self.alphabet_b = contents[7]
        for i in range(8, len(contents)):
            temp = contents[i].split()
            if len(temp) > 0:
                self.match_matrix.set_score(temp[2], temp[3], float(temp[4]))
        
class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """

    # Function to initialize the Align class
    # Inputs:
    #       None
    # Returns:
    #       None
    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters() 

    # Function that drives the sequence alignment logic
    # Inputs:
    #       None
    # Returns:
    #       None
    def align(self):
        """
        Main method for running alignment.
        """
        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # perform a traceback and write the output to an output file
        max_val, max_loc = self.find_traceback_start()
        print(max_val)
        all_alignments = []
        for i in max_loc:
            self.traceback("M", i[0], i[1], "", "", "", all_alignments)
        self.write_output(self.output_file, max_val, all_alignments)

    # Function that sweeps through the matrices and updates their states
    # Inputs:
    #       None
    # Returns:
    #       None
    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        """
        self.m_matrix = ScoreMatrix("M", self.align_params.num_rows, self.align_params.num_cols)
        self.ix_matrix = ScoreMatrix("Ix", self.align_params.num_rows, self.align_params.num_cols)
        self.iy_matrix = ScoreMatrix("Iy", self.align_params.num_rows, self.align_params.num_cols)

        for i in range(1, self.align_params.num_rows):
            for j in range(1, self.align_params.num_cols):
                self.update(i, j)

    # Function that calls the update function for all three matrices
    # Inputs:
    #       row: current row to be updated
    #       col: current col to be updated
    # Returns:
    #       None
    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.

        Input:
           row = the row index to update
           col = the column index to update
        """
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    # Function that updates the M matrix at a given row and column
    # Inputs:
    #       row: current row to be updated
    #       col: current col to be updated
    # Returns:
    #       None
    def update_m(self, row, col):
        match_ij = self.align_params.match_matrix.get_score(self.align_params.seq_a[row - 1], self.align_params.seq_b[col - 1])
        no_gap = self.m_matrix.get_score(row - 1, col - 1) + match_ij
        gap_b = self.ix_matrix.get_score(row - 1, col - 1) + match_ij
        gap_a = self.iy_matrix.get_score(row - 1, col - 1) + match_ij
        score = max(no_gap, gap_b, gap_a)
        if self.align_params.global_alignment:
            self.m_matrix.set_score(row, col, score) 
        else:
            self.m_matrix.set_score(row, col, max(0, score))

        if fuzzy_equals(score, no_gap):
            self.m_matrix.set_pointers(row, col, ("M", row - 1, col - 1))
        if fuzzy_equals(score, gap_b):
            self.m_matrix.set_pointers(row, col, ("Ix", row - 1, col - 1))
        if fuzzy_equals(score, gap_a):
            self.m_matrix.set_pointers(row, col, ("Iy", row - 1, col - 1))

    # Function that updates the Ix matrix at a given row and column
    # Inputs:
    #       row: current row to be updated
    #       col: current col to be updated
    # Returns:
    #       None
    def update_ix(self, row, col):
        first_gap = self.m_matrix.get_score(row - 1, col) - self.align_params.dy
        additional_gap = self.ix_matrix.get_score(row - 1, col) - self.align_params.ey
        score = max(first_gap, additional_gap)
        if self.align_params.global_alignment:
            self.ix_matrix.set_score(row, col, score)
        else:
            self.ix_matrix.set_score(row, col, max(0, score))

        if fuzzy_equals(score, first_gap):
            self.ix_matrix.set_pointers(row, col, ("M", row - 1, col))
        if fuzzy_equals(score, additional_gap):
            self.ix_matrix.set_pointers(row, col, ("Ix", row - 1, col))

    # Function that updates the Iy matrix at a given row and column
    # Inputs:
    #       row: current row to be updated
    #       col: current col to be updated
    # Returns:
    #       None
    def update_iy(self, row, col):
        first_gap = self.m_matrix.get_score(row, col - 1) - self.align_params.dx
        additional_gap = self.iy_matrix.get_score(row, col - 1) - self.align_params.ex
        score = max(first_gap, additional_gap)
        if self.align_params.global_alignment:
            self.iy_matrix.set_score(row, col, score)
        else:
            self.iy_matrix.set_score(row, col, max(0, score))

        if fuzzy_equals(score, first_gap):
            self.iy_matrix.set_pointers(row, col, ("M", row, col - 1))
        if fuzzy_equals(score, additional_gap):
            self.iy_matrix.set_pointers(row, col, ("Iy", row, col - 1))

    # Function that finds the best alignment score and where to start the traceback
    # Inputs:
    #       None
    # Returns:
    #       the best alignment score and a set of coordinates to start the traceback from
    def find_traceback_start(self):
        """
        Finds the location to start the traceback..
        Think carefully about how to set this up for local 

        Returns:
            (max_val, max_loc) where max_val is the best score
            max_loc is a set() containing tuples with the (i,j) location(s) to start the traceback
             (ex. [(1,2), (3,4)])
        """
        max_val = 0.0
        max_loc = set()
        if self.align_params.global_alignment:
            for i in range(self.align_params.num_cols):
                if self.m_matrix.get_score(self.align_params.num_rows - 1, i) > max_val:
                    max_val = self.m_matrix.get_score(self.align_params.num_rows - 1, i)
                    max_loc = set([(self.align_params.num_rows - 1, i)])
                elif fuzzy_equals(self.m_matrix.get_score(self.align_params.num_rows - 1, i), max_val):
                    max_loc.add((self.align_params.num_rows - 1, i))
            for i in range(self.align_params.num_rows):
                if self.m_matrix.get_score(i, self.align_params.num_cols - 1) > max_val:
                    max_val = self.m_matrix.get_score(i, self.align_params.num_cols - 1)
                    max_loc = set([(i, self.align_params.num_cols - 1)])
                elif fuzzy_equals(self.m_matrix.get_score(i, self.align_params.num_cols - 1), max_val):
                    max_loc.add((i, self.align_params.num_cols - 1))
        else:
            for i in range(self.align_params.num_rows):
                for j in range(self.align_params.num_cols):
                    if self.m_matrix.get_score(i, j) > max_val:
                        max_val = self.m_matrix.get_score(i, j)
                        max_loc = set([(i, j)])
                    elif fuzzy_equals(self.m_matrix.get_score(i, j), max_val):
                        max_loc.add((i, j))
        return max_val, max_loc

    # Function that recursively reconstructs all alignments that produce the best score
    # Inputs:
    #       matrix: current matrix
    #       row: current row
    #       col: current column
    #       seq_a: sequence a in progress
    #       seq_b: sequence b in progress
    #       path: path of traceback in progress (used for debugging)
    #       all_alignments: list holding all completed alignments
    # Returns:
    #       None
    def traceback(self, matrix, row, col, seq_a, seq_b, path, all_alignments):
        if self.align_params.global_alignment:
            if row == 0 or col == 0:
                if not (seq_a, seq_b) in all_alignments:
                    all_alignments.append((seq_a, seq_b))
                    return
        else:
            if matrix == "M" and self.m_matrix.get_score(row, col) == 0:
                if not (seq_a, seq_b) in all_alignments:
                    all_alignments.append((seq_a, seq_b))
                    return
            elif matrix == "Ix" and self.ix_matrix.get_score(row, col) == 0 or matrix == "Iy" and self.iy_matrix.get_score(row, col) == 0:
                return
                
        if matrix == "M":
            new_a = self.align_params.seq_a[row - 1] + seq_a
            new_b = self.align_params.seq_b[col - 1] + seq_b
        elif matrix == "Ix":
            new_a = self.align_params.seq_a[row - 1] + seq_a
            new_b = "_" + seq_b
        elif matrix == "Iy":
            new_a = "_" + seq_a
            new_b = self.align_params.seq_b[col - 1] + seq_b
            
        if matrix == "M":
            pointers = self.m_matrix.get_pointers(row, col)
        elif matrix == "Ix":
            pointers = self.ix_matrix.get_pointers(row, col)
        elif matrix == "Iy":
            pointers = self.iy_matrix.get_pointers(row, col)

        for i in pointers:
            self.traceback(i[0], i[1], i[2], new_a, new_b, path + "->" + str(i), all_alignments)

    # Function that writes the best alignment score and all alignments to a file
    # Inputs:
    #       output_file: name of file to write outputs to
    #       max_val: best alignment score
    #       all_alignments: list holding all completed alignments
    # Returns:
    #       None
    def write_output(self, output_file, max_val, all_alignments):
        with open(output_file, "w") as f:
            f.write(str(round(max_val, 1)))
            f.write("\n\n")
            for i in all_alignments:
                f.write(i[0] + "\n")
                f.write(i[1] + "\n\n")

# Function that reads from the command line and instantiates the Align class
# Inputs:
#       None
# Returns:
#       None
def main():
    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()

if __name__=="__main__":
    main()