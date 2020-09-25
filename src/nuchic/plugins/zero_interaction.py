import logger as logger
from nuchic.pyinteractions import PyInteraction
import interactions as interactions
import vectors


class ZeroInteraction(PyInteraction):
    """ Interaction model that always returns a constant value. """
    _registered = True
    name = 'ZeroInteraction'

    def __init__(self, *args):
        # Remove the name argument
        args = args[1:]

        interactions.Interactions.__init__(self)
        logger.info(f"ZeroInteraction: xsec = 0")
        super().__init__()

    def CrossSection(self, part1, part2):  # pylint: disable=unused-argument
        """ Return the cross-section value. """
        return 0

    def IsRegistered(self):  # pylint: disable=invalid-name
        """ Check if the interaction has been registered in the factory. """
        return self._registered

    def MakeMomentum(self, same_pid, p1_cm, pcm, rans):  # pylint: disable=unused-argument
        """ Make the momentum for the outgoing nucleons. """
        return vectors.Vector3()
