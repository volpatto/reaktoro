# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2015 Allan Leal
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import yaml
from reaktoro import *


# A type used to define operations that interpret a Reaktoro script file.
class Interpreter:

    # Construct a default an Interpreter instance
    def __init__():
        # The database used to initialize the chemical system
        self.database = None

        # The chemical editor used to define the chemical system
        self.editor = None

        # The chemical system used for the calculations
        self.system = None

        # The chemical reactions controlled by kinetics
        self.reactions = None

        # The map of chemical states produced during the calculation
        self.states = {}

        # The defined mineral reactions
        self.mineral_reactions = []

        # The list of all compound names found in the script file.
        self.compounds = []

        # The list of all elements that compose the compounds.
        self.elements = []

    # Initialize essential members of the interpreter state.
    # This method is used to initialize a ChemicalEditor instance by firstly
    # identifying all compound and species names in the input file. After this,
    # an aqueous phase is created with all possible species found in the database
    # that contains the elements composing the list of found compounds. A gaseous
    # phase is created only if names of gaseous species present in the database are
    # found in the input script. Pure mineral phases are created for each mineral name
    # found in the input script file that is also present in the database.
    def init(root):
        # Initialize the database
        self.database = Database('supcrt98')

        # Initialize the list of compounds found in the script file
        self.compounds = collectCompounds(root)

        # Initialize the list of elements that compose the compounds
        self.elements = identifyElements(self.compounds, self.database)

        # Initialize the automatic chemical system
        initChemicalEditor(root)

    # Execute a Reaktoro input script as string.
    def execute(sript):
        # Preprocess the input script so that it conforms with YAML rules.
        str = preprocess(sript)

        # The root node of the yaml script
        root = yaml.load(str)

        # Initialize some essential members of the interpreter state
        init(root)

        # The map of process functions (from keyword to respective process function)
        processfn_map = {
            'database': processDatabaseNode,
            'aqueousphase': processAqueousPhaseNode,
            'gaseousphase': processGaseousPhaseNode,
            'mineralphase': processMineralPhaseNode,
            'minerals': processMineralsNode,
            'mineralreaction': processMineralReactionNode,
            'equilibrium': processEquilibriumNode,
            'equilibriumproblem': processEquilibriumNode,
            'equilibriumpath': processEquilibriumPathNode,
            'kinetics': processKineticPathNode,
            'kineticpath': processKineticPathNode,
            'phreeqc': processPhreeqcNode,
        }

        # Keyword triggers that start the initialization of chemical system
        triggers = [
            'equilibrium',
            'equilibriumproblem',
            'speciation',
            'speciationproblem',
        ]

        # For every child node in the root node...
        for node in root:
            # The keyword of the current yaml node
            key = lowercase(keyword(node))

            # Check if current key is a chemical system initializer
            if key in triggers and self.system is None:
                self.system = ChemicalSystem(self.editor)

            # Find an entry in the process function map with that key
            processfn = processfn_map.get(key)

            # Assert that this entry was found
            assert processfn is not None, 'Could not parse `' + node + '`. ' \
                'Expecting a valid keyword. Did you misspelled it?'

            # Process the current child node and update the interpreter state
            processfn(node)

    def processDatabaseNode(node, identifier):
        self.database = Database(str(valnode(node)))
        self.editor = ChemicalEditor(self.database)
        self.elements = identifyElements(self.compounds, self.database)

    def processAqueousPhaseNode(node, identifier):
        phasename = identifier
        assert phasename is None, 'Currently, the name of the aqueous phase ' \
            'cannot be changed from `Aqueous` to %s' % identifier
        compounds = node
        if compounds == 'auto':
            self.editor.addAqueousPhaseWithElements(self.elements)
        else self.editor.addAqueousPhase(compounds)

    def processGaseousPhaseNode(node, identifier):
        phasename = identifier
        assert phasename is None, 'Currently, the name of the gaseous phase ' \
            'cannot be changed from `Gaseous` to %s' % identifier

        assert phasename, 'Could not set the gaseous phase with `' + node + '`.' \
            'The name of the gaseous phase is `Gaseous` and cannot be changed.'
        compounds = node
        if compounds == 'auto':
            self.editor.addGaseousPhaseWithElements(self.elements)
        else self.editor.addGaseousPhase(compounds)

    def processMineralPhaseNode(node, identifier):
        phasename = identifier
        mineralname = node
        phase = self.editor.addMineralPhase(mineralname)
        if phasename:
            phase.setName(phasename)

    def processMineralsNode(node, identifier):
        minerals = node
        if minerals == 'auto':
            for mineral in self.database.mineralSpeciesWithElements(self.elements):
                self.editor.addPhase(MineralPhase(mineral))
        else for mineral in minerals:
            self.editor.addMineralPhase(mineral)

    def processChemicalModelNode(node, identifier):
        pass

    def processMineralReactionNode(node, identifier):
        # Convert the yaml node into a keyword
        kwd::MineralReaction keyword; node >> keyword

        # Convert the keyword into a MineralReaction instance
        reaction = MineralReaction()
        initializeMineralReaction(reaction, keyword)

        # Add the new mineral reaction to the list
        self.mineral_reactions.append(reaction)

    def processEquilibriumNode(node, identifier):
        # Convert the yaml node into a keyword
        kwd::EquilibriumProblem keyword; node >> keyword

        # Initialize the equilibrium problem using the just initialized chemical system
        EquilibriumProblem problem(self.system)
        initializeEquilibriumProblem(problem, keyword)

        # Initialize the chemical state
        state = EquilibriumState(self.system)

        # Initialize the amounts of the inert species
        for s in keyword.inert_species:
            state.setSpeciesAmount(s.entity, s.value, s.units)

        # Perform the equilibrium calculation
        equilibrate(state, problem)

        # Output the resulting chemical state
        state.output(keyword.stateid + '.dat')

        # Store the resulting chemical state with given state ID
        self.states[keyword.stateid] = state

    def processEquilibriumPathNode(node, identifier):
        # Convert the yaml node into a kinetic keyword
        kwd::EquilibriumPath keyword; node >> keyword

        # Initialize the kinetic path instance
        path = EquilibriumPath(self.system)

        # Set the inert species in the equilibrium path calculation
        partition = Partition(system)
        partition.setInertSpecies(keyword.inert_species)
        path.setPartition(partition)

        # Set the plots in the kinetic path problem
        for plot_options in keyword.plots:
            plot = path.plot()
            initializeChemicalPlot(plot, plot_options)


        initializeEquilibriumPath(path, keyword)

        # Alias to the initial and final states
        statei = self.states[keyword.initial_state]
        statef = self.states[keyword.final_state]

        # Solve the equilibrium path problem
        path.solve(statei, statef)

    def processKineticPathNode(node, identifier):
        # Convert the yaml node into a kinetic keyword
        kwd::KineticPath keyword; node >> keyword

        # Initialize the ReactionSystem instance
        for reaction in self.mineral_reactions:
            self.editor.addMineralReaction(reaction)

        # Set the reactions of the chemical system
        self.reactions = self.editor

        # Initialize the kinetic path instance
        path = KineticPath(self.reactions)
        initializeKineticPath(path, keyword)

        # Alias to the initial condition state
        auto& state = self.states[keyword.initial_condition]

        # The duration time and its units
        duration = keyword.duration.value
        units = keyword.duration.units

        # Solve the kinetics problem
        path.solve(state, 0, duration, units)

        # Output the final chemical state with given state id
        state.output(keyword.stateid + '.dat')

        # Store the final chemical state into the map of states
        self.states[keyword.stateid] = state

    def processPhreeqcNode(node, identifier):
        # Convert the yaml node into a kinetic keyword
        kwd::PhreeqcKeyword keyword; node >> keyword

        # Initialize a Phreeqc instance with given database and input script
        phreeqc = Phreeqc()
        phreeqc.load(keyword.database)
        phreeqc.execute(keyword.input, keyword.output)

        # Initialize the chemical system
        self.system = ChemicalSystem(phreeqc)

        # Initialize the chemical state with the current state of phreeqc
        state = phreeqc.state(self.system)

        # Output the final chemical state with given state id
        state.output(keyword.stateid + '.dat')

        # Store the final chemical state into the map of states
        self.states[keyword.stateid] = state


#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

// Reaktoro includes
#include <unsupported/cpp-interpreter/Interpreter.hpp>
#include <unsupported/cpp-interpreter/Processors.hpp>
#include <unsupported/cpp-interpreter/Utils.hpp>

namespace Reaktoro {

struct Interpreter::Impl
{
    /// The state of the interpreter
    InterpreterState istate

    /// Construct a default Impl instance
    Impl()
    {}



    // Initialize the default state of chemical editor based on the compounds found in the script
    auto initChemicalEditor(const Node& root) -> void
    {
        // Auxiliary references
        const auto& database = self.database
        const auto& elements = self.elements

        // Collect all aqueous species that can be formed out of the elements found in the script
        auto aqueous_species = names(database.aqueousSpeciesWithElements(elements))

        // Determine if there are gaseous and mineral species among the found compounds in the script
        auto mineral_species = filterMineralSpecies(self.compounds, self.database)
        auto gaseous_species = filterGaseousSpecies(self.compounds, self.database)

        // Check if there is any Speciation, in which case all possible gases and minerals are considered
        if(hasSpeciation(root))
        {
            gaseous_species = names(database.gaseousSpeciesWithElements(elements))
            mineral_species = names(database.mineralSpeciesWithElements(elements))
        }

        // Add the aqueous phase using the compounds as the initializer
        self.editor.addAqueousPhaseWithSpecies(aqueous_species)

        // Add a gaseous phase if there are gaseous species
        if(gaseous_species.size())
            self.editor.addGaseousPhaseWithSpecies(gaseous_species)

        // Add a mineral phase for each mineral species
        for(auto x : mineral_species)
            self.editor.addMineralPhaseWithSpecies({x})
    }

