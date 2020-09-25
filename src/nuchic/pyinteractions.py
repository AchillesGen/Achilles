import interactions

class PyInteraction(interactions.Interactions):
    _registered = True

    class Unknown(Exception):
        """ Error to raise for unknown RunMode. """

    @classmethod
    def get_all_subclasses(cls):
        """ Get all the subclasses recursively of this class. """
        for subclass in cls.__subclasses__():
            yield from subclass.get_all_subclasses()
            yield subclass

    @classmethod
    def list_all(cls):
        """ List all the subclasses. """
        print('Python Interactions:')
        print('--------------------')
        for subclass in cls.get_all_subclasses():
            print(f' - {subclass.name}')
        print()

    def __new__(cls, name, *args, **kwargs):
        del args
        del kwargs
        for subclass in cls.get_all_subclasses():
            if subclass.name == name:
                # Using "object" base class methods avoids recursion here.
                return interactions.Interactions.__new__(subclass)
        raise PyInteraction.Unknown('Interaction "{}" is not known.'.format(name))

    def __init__(self):
        """ General run mode initialization. """
        instance_name = self.__class__.__name__

    def CrossSection(self, part1, part2):
        """ Calculate the cross section. """
        raise NotImplementedError

    def MakeMomentum(self, events):
        """ Make the final state momentums. """
        raise NotImplementedError
