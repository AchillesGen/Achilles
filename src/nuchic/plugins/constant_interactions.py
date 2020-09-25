import logger as logger
from nuchic.pyinteractions import PyInteraction
import interactions as interactions
import vectors


class ConstantInteraction(PyInteraction):
    """ Interaction model that always returns a constant value. """
    _registered = True
    name = 'ConstantInteraction'

    def __init__(self, *args):
        # Remove the name argument
        args = args[1:]

        interactions.Interactions.__init__(self)
        self.xsec = args[0]
        logger.info(f"ConstantInteraction: xsec = {self.xsec}")
        super().__init__()

    def CrossSection(self, part1, part2):  # pylint: disable=unused-argument
        """ Return the cross-section value. """
        return self.xsec

    def IsRegistered(self):  # pylint: disable=invalid-name
        """ Check if the interaction has been registered in the factory. """
        return self._registered

    def MakeMomentum(self, same_pid, p1_cm, pcm, rans):  # pylint: disable=unused-argument
        """ Make the momentum for the outgoing nucleons. """
        return vectors.Vector3()
