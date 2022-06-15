# This file was automatically created by FeynRules 2.3.36
# Mathematica version: 12.0.0 for Microsoft Windows (64-bit) (April 6, 2019)
# Date: Tue 10 Aug 2021 11:00:24


from .object_library import all_orders, CouplingOrder


QCD = CouplingOrder(name = 'QCD',
                    expansion_order = 99,
                    hierarchy = 1)

QED = CouplingOrder(name = 'QED',
                    expansion_order = 99,
                    hierarchy = 2)

NP = CouplingOrder(name = 'NP',
                   expansion_order = 99,
                   hierarchy = 1)

NucPhys = CouplingOrder(name = 'NucPhys',
                        expansion_order = 1,
                        hierarchy = 1)
