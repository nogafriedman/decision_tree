import itertools


class Node:
    def __init__(self, data, positive_child=None, negative_child=None):
        self.data = data  # question asked, final decision if leaf
        self.positive_child = positive_child  # positive answer, none if leaf
        self.negative_child = negative_child  # negative answer, none if leaf


class Record:
    def __init__(self, illness, symptoms):
        self.illness = illness  # str - name of illness
        self.symptoms = symptoms  # list of strings - symptoms


def parse_data(filepath):
    with open(filepath) as data_file:
        records = []
        for line in data_file:
            words = line.strip().split()
            records.append(Record(words[0], words[1:]))
        return records


class Diagnoser:
    def __init__(self, root):
        self.root = root

    def diagnose(self, symptoms):  # 1
        """
		Recursive version of diagnose.
		:param symptoms: list of symptoms as strings
		:return: str - diagnosis based on the symptoms
		"""
        node = self.root
        return self._diagnose_helper(node, symptoms)

    def _diagnose_helper(self, node, symptoms):  # recursive func for diagnose
        if node.positive_child is None or node.negative_child is None:
            diagnosis = node.data
            return diagnosis
        if node.data in symptoms:
            return self._diagnose_helper(node.positive_child, symptoms)
        else:
            return self._diagnose_helper(node.negative_child, symptoms)

    def calculate_success_rate(self, records):  # 2
        if len(records) == 0:
            raise ValueError("Records is an empty list!")
        success = 0
        for rec in records:
            diagnosis = self.diagnose(rec.symptoms)
            if diagnosis == rec.illness:
                success += 1
        res = success / len(records)
        return res

    # def diagnose(self, symptoms):  # 1
    # 	"""
    # 	Loop version of diagnose.
    # 	:param symptoms: list of symptoms as strings
    # 	:return: str - diagnosis based on the symptoms
    # 	"""
    # 	node = self.root
    # 	# return self.diagnose_helper(node, symptoms)
    # 	while node.positive_child is not None:
    # 		if node.data in symptoms:
    # 			node = node.positive_child
    # 		else:
    # 			node = node.negative_child
    # 	return node.data

    def all_illnesses(self):  # 3
        """
		Goes through the tree and returns a list with all of the illnesses
		listed in it (aka - the leaves of the tree), organized by number of
		appearances in the tree, from most common to least.
		:return: list
		"""
        illnesses = dict()
        node = self.root
        # go to recursive function:
        self._illnesses_helper(node, illnesses)
        # sort the returned dictionary by most common to least:
        sorted_list = sorted(illnesses.items(), key=lambda x: x[1],
                             reverse=True)
        illnesses_list = []
        # add illnesses from dictionary to list after sorted:
        for illness in sorted_list:
            illnesses_list.append(illness[0])
        return illnesses_list

    def _illnesses_helper(self, node, illnesses):
        """
		Recursively goes through the tree and adds the illness to the
		dictionary only if it's a leaf, or adds +1 to it's counter if it's
		already in the dictionary.
		:param node: current node
		:param illnesses: dictionary with illness name (key), number of
						  appearances on the tree (value)
		:return: updated dictionary
		"""
        if node.positive_child is None:  # stopping condition: node is a leaf
            if node.data is not None:
                if node.data in illnesses:
                    illnesses[node.data] += 1  # if illness already in dictionary
                else:
                    illnesses[node.data] = 1  # if illness is not yet in dictionary
            return
        self._illnesses_helper(node.positive_child, illnesses)
        self._illnesses_helper(node.negative_child, illnesses)
        return

    def paths_to_illness(self, illness):  # 4
        """
		Receives an illness name and returns all the paths leading to it in
		the tree.
		:param illness: str
		:return: list of lists, each sublist contains True/False according
		to the turns taken in the tree
		"""
        node = self.root
        path_list = []
        if node.data is None:
            return path_list
        self._path_helper(node, illness, path_list, list())
        return path_list

    def _path_helper(self, node, illness, path_list, cur_path):
        """
		Recursively goes through the tree and adds True/False o the list
		according to the turns taken. Adds path to list only if the final
		leaf was the illness it looked for.
		:param node: node object
		:param illness: str
		:param path_list: list (to contain all paths)
		:param cur_path: list (to append to path_list)
		:return: None
		"""
        if node.positive_child is None:  # stopping condition
            if node.data == illness:  # add only if leaf is illness
                path_list.append(cur_path)
            return
        # try going right:
        pos_path = cur_path.copy()
        pos_path.append(True)
        self._path_helper(node.positive_child, illness, path_list, pos_path)
        # try going left:
        neg_path = cur_path.copy()
        neg_path.append(False)
        self._path_helper(node.negative_child, illness, path_list, neg_path)
        return

    def minimize(self, remove_empty=False):  # 7
        """
        If function is called with parameter remove_empty as False:
        Removes all redundant nodes from the tree - a redundant node is a
        node whose paths reached both the same diagnosis.
        If function is called with parameter remove_empty as True:
        In addition to the changes made for False, also removes leaves
        whose data is None, and replaces their parent node with the second
        child.
        :param remove_empty: boolean - True or False
        :return: None
        """
        # check if root is leaf:
        if self.is_leaf(self.root):
            return
        self._minimize_helper(remove_empty, self.root)
        return

    def _minimize_helper(self, remove_empty, node, parent_node=None):
        """
        Recursively checks for redundant nodes in the tree and removes them.
        :param remove_empty: boolean - True/False
        :param node: node object (starts as the root)
        :param parent_node: parent node for the current node (default = None)
        :return: None
        """
        if self.is_leaf(node):
            # we reached the leaf, recurse back:
            return
        if node.positive_child.data is not None:
            self._minimize_helper(remove_empty, node.positive_child, node)
        if node.negative_child.data is not None:
            self._minimize_helper(remove_empty, node.negative_child, node)
        if remove_empty:
            self._remove_true(node, parent_node)
        is_redundant, leaves = self._check_redundancy(node)
        if is_redundant:
            self._remove_false(node, parent_node, leaves)

    def _check_redundancy(self, root):
        """
        Gets the leaves the converging paths and compares them. If the
        leaves are the same then their root is redundant.
        :param root: node object
        :return: True if node was found redundant, False if not
        """
        pos_leaves = list()
        neg_leaves = list()
        if root.positive_child.data is not None:
            self._get_leaves(root.positive_child, pos_leaves)
        if root.negative_child.data is not None:
            self._get_leaves(root.negative_child, neg_leaves)
        leaves = [pos_leaves, neg_leaves]
        if pos_leaves == neg_leaves:
            return True, leaves
        return False, None

    def _remove_false(self, node, parent_node, leaves):
        if len(leaves[0]) > 1:
            if parent_node is None:  # node is the root of the tree
                self.root = node.positive_child
                return
            else:
                new_node = node.positive_child
                if parent_node.positive_child is node:
                    parent_node.positive_child = new_node
                if parent_node.negative_child is node:
                    parent_node.negative_child = new_node
                    return
        else:
            new_node = node.positive_child
            if parent_node is None:  # node is the root of the tree
                self.root = new_node
                return
            if parent_node.positive_child is node:
                parent_node.positive_child = new_node
            if parent_node.negative_child is node:
                parent_node.negative_child = new_node
            return

    def _remove_true(self, node, parent_node):
        if node.positive_child.data is None:
            new_child = node.negative_child
            if parent_node is None:
                self.root = new_child
                return
            else:
                if parent_node.positive_child is node:
                    parent_node.positive_child = new_child
                if parent_node.negative_child is node:
                    parent_node.negative_child = new_child
        if node.negative_child.data is None:
            new_child = node.positive_child
            if parent_node is None:
                self.root = new_child
                return
            else:
                if parent_node.positive_child is node:
                    parent_node.positive_child = new_child
                if parent_node.negative_child is node:
                    parent_node.negative_child = new_child

    def _get_leaves(self, node, leaves):
        """
        Goes through the tree and get all the leaves from positive to negative
        :param node: node object
        :param leaves: list of strings - string is the leaf's data
        :return:
        """
        if node.positive_child is not None:
            self._get_leaves(node.positive_child, leaves)
        if node.negative_child is not None:
            self._get_leaves(node.negative_child, leaves)
        if self.is_leaf(node):
            leaves.append(node.data)

    def is_leaf(self, node):
        """
        Checks if the node is a leaf.
        :param node: node object
        :return: True if leaf, False if not
        """
        if node.positive_child is None and node.negative_child is None:
            return True
        return False


def build_tree(records, symptoms):  # 5
    """
	Builds a tree that asks about the symptoms received as parameter
	:param records: list of objects of type Record
	:param symptoms: list of strings
	:return: object of class Diagnoser
	"""
    # deal with exceptions:
    _check_type(records, symptoms)
    if len(symptoms) == 0:
        illness = _find_diagnosis(records, symptoms, [])
        diagnoser = Diagnoser(Node(illness))
        return diagnoser
    # create a root for the tree:
    root = Node(symptoms[0])
    # create a new diagnoser object with root:
    diagnoser = Diagnoser(root)
    _build_helper(records, symptoms, root, list())
    return diagnoser


def _build_helper(records, symptoms, root, path_bool, ind=0):
    if ind == len(symptoms) - 1:
        # this means current node is the last symptom on the list,
        # create two leaf children and calculate each one's matching illness:
        path_copy = path_bool.copy()
        path_bool.append(True)
        pos_illness = _find_diagnosis(records, symptoms, path_bool)
        pos_leaf = Node(pos_illness, None, None)
        path_copy.append(False)
        neg_illness = _find_diagnosis(records, symptoms, path_copy)
        neg_leaf = Node(neg_illness, None, None)
        root.positive_child = pos_leaf
        root.negative_child = neg_leaf
        return
    # update data for current root object:
    pos_node = Node(symptoms[ind + 1])
    neg_node = Node(symptoms[ind + 1])
    root.positive_child = pos_node
    root.negative_child = neg_node
    # positive route:
    pos_path = path_bool.copy()
    pos_path.append(True)
    _build_helper(records, symptoms, pos_node, pos_path, ind + 1)
    # negative route:
    neg_path = path_bool.copy()
    neg_path.append(False)
    _build_helper(records, symptoms, neg_node, neg_path, ind + 1)


def _find_diagnosis(records, symptoms, path_bool):
    """
    Sorts the bool list created in the recursion into a list of True
    symptoms (that were answered 'yes' to) and False symptoms (that were
    answered 'no' to). Then goes through the records and checks if any
    record matches the symptoms - if so, adds the diagnosis to a dictionary.
    Finally, return the most common diagnosis to the build_tree function.
    :param records: list of objects of class Record
    :param symptoms: list of strings
    :param path_bool: list of booleans - True/False
    :return: string - illness
    """
    # create a list for symptoms that were True and for ones that were False:
    lst = _split_bools(symptoms, path_bool)
    symptoms_true = lst[0]
    symptoms_false = lst[1]

    # compare lists to records to find the matching diagnosis:
    possible_diagnosis = {}
    for record in records:
        if _check_true_symptoms(record.symptoms, symptoms_true) and \
                _check_false_symptoms(record.symptoms, symptoms_false):
            if record.illness in possible_diagnosis:
                possible_diagnosis[record.illness] += 1
            else:
                possible_diagnosis[record.illness] = 1
    # sort the returned dictionary by most common illness to least:
    sorted_list = sorted(possible_diagnosis.items(), key=lambda x: x[1],
                         reverse=True)
    if not sorted_list:
        return None
    return sorted_list[0][0]


def _split_bools(symptoms, path_bool):
    """
	Receives a list of bools and splits them into a True list and a False list.
	:param symptoms: list of symptoms
	:param path_bool: list with element True/False
	:return: list of two lists - index 1 is True and index 2 is False.
	"""
    symptoms_true = []
    symptoms_false = []
    for symptom_ind in range(len(symptoms)):
        if path_bool[symptom_ind]:
            symptoms_true.append(symptoms[symptom_ind])
        else:
            symptoms_false.append(symptoms[symptom_ind])
    return [symptoms_true, symptoms_false]


def _check_true_symptoms(record, symptoms_true):
    """
	Compares the received symptoms list to the record list and returns True
	only if every symptom in the symptoms list really appears in the record
	list.
	:param record: list of symptoms
	:param symptoms_true: list of symptoms
	:return: True if the lists matched, False if not
	"""
    match = True
    for symptom in symptoms_true:
        if symptom not in record:
            match = False
            break  # it's not this record! try the next one
        else:
            match = True
    return match


def _check_false_symptoms(record, symptoms_false):
    """
	Compares the received symptoms list to the record list and returns True
	only if every symptom in the symptoms list DIDN'T appear in the record
	list.
	:param record: list of symptoms
	:param symptoms_false: list of symptoms
	:return: True if the lists matched, False if not
	"""
    match = True
    for symptom in symptoms_false:
        if symptom in record:
            match = False
            break  # it's not this record! try the next one
        else:
            match = True
    return match


def _check_type(records, symptoms):
    """
    Checks the parameters are valid and if not, raises an exception.
    :param records: list of objects of class Record
    :param symptoms: list of strings
    :return: None
    """
    for record in records:
        if not isinstance(record, Record):
            raise TypeError("All elements of records parameter should be "
                            "record objects!")
    for symptom in symptoms:
        if not isinstance(symptom, str):
            raise TypeError("All elements of symptoms parameter should be "
                            "strings!")


def optimal_tree(records, symptoms, depth):  # 6
    """
    Creates all subset variations of the symptoms list with len = depth,
    then build a tree for it using build_tree.
    Finally, calculate the success rate of that tree using
    calculate_success_rate function (#2 function in this file), and returns
    the tree (Diagnoser object) of the highest success rate.
    :param records: list of objects of class Record
    :param symptoms: list of strings
    :param depth: int
    :return: object of class Diagnoser
    """
    check_exceptions_optimal(records, symptoms, depth)  # catch exceptions
    if depth == 0:
        optimal = build_tree(records, [])
        return optimal
    if len(records) == 0:
        trees_list = []
        for subset in itertools.combinations(symptoms, depth):
            diagnoser = build_tree(records, subset)
            trees_list.append(diagnoser)
        return trees_list[0]
    else:
        trees_list = []
        for subset in itertools.combinations(symptoms, depth):
            diagnoser = build_tree(records, subset)
            success_rate = diagnoser.calculate_success_rate(records)
            trees_list.append((diagnoser, success_rate))
        trees_list.sort(key=lambda x: x[1], reverse=True)
        optimal = trees_list[0][0]
    return optimal


def check_exceptions_optimal(records, symptoms, depth):
    if depth > len(symptoms):
        raise ValueError("Depth cannot be longer than the symptoms list!")
    if depth < 0:
        raise ValueError("Depth has to be a non-negative int!")
    if len(symptoms) != len(set(symptoms)):
        raise ValueError("Each symptom must appear on the list only once!")
    for record in records:
        if not isinstance(record, Record):
            raise TypeError("All record elements should be of type Record!")
    for symptom in symptoms:
        if not isinstance(symptom, str):
            raise TypeError("All symptoms should be of type str!")
