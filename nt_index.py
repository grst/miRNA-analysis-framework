import warnings

class QuickNTIndex:
    """ create an index of fixed length
    over an input string so that search is possible
    in constant time
    """
    def __init__(self, haystack, iLen=6):
        self.haystack = haystack
        self.iLen = iLen
        self._make_index()

    def _make_index(self):
        self.indexmap = {}
        for i in range(len(self.haystack) - self.iLen + 1):
            tmpIndex = self.haystack[i:i+self.iLen]
            if(self.indexmap.get(tmpIndex) == None):
                self.indexmap[tmpIndex] = []
            self.indexmap[tmpIndex].append(i)

    def index(self, needle):
        if len(needle) < self.iLen:
            warnings.warn("needle shorter than index, string will not be found: " + needle)
        indexes = self.indexmap.get(needle[:self.iLen])
        if indexes == None: return -1
        realIndexes = []
        for i in indexes:
            if(self.haystack.find(needle, i, i+len(needle)) >= 0):
                return i
        return -1

    def all_indexes(self, needle):
        if(len(needle) < self.iLen):
            warnings.warn("needle shorter than index, string will not be found: " + needle)
        indexes = self.indexmap.get(needle[:self.iLen])
        if indexes == None: return []
        realIndexes = []
        for i in indexes:
            if(self.haystack.find(needle, i, i+len(needle)) >= 0):
                realIndexes.append(i)
        return realIndexes
