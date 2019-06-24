__author__ = "Daniel Winklehner"


class MyColors(object):

    def __init__(self):
        """
        Constructor
        """
        self.colors = ['#4B82B8',
                       '#B8474D',
                       '#95BB58',
                       '#234B7C',
                       '#8060A9',
                       '#53A2CB',
                       '#FC943B']

    def __getitem__(self, item):

        return self.colors[int(item % 7)]
