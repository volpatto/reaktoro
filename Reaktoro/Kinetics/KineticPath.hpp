// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <memory>
#include <string>
#include <vector>

namespace Reaktoro {

// Forward declarations
class ChemicalOutput;
class ChemicalPlot;
class ChemicalState;
class Partition;
class ReactionSystem;
struct KineticOptions;

/// A class that conveniently solves kinetic path calculations.
class KineticPath
{
public:
    /// Construct a KineticPath instance.
    explicit KineticPath(const ReactionSystem& reactions);

    /// Construct a copy of a KineticPath instance.
    KineticPath(const KineticPath& other);

    /// Destroy the KineticPath instance.
    virtual ~KineticPath();

    /// Assign a KineticPath instance to this instance.
    auto operator=(KineticPath other) -> KineticPath&;

    /// Set the options for the chemical kinetics calculation.
    auto setOptions(const KineticOptions& options) -> void;

    /// Set the partition of the chemical system.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    auto setPartition(const Partition& partition) -> void;

    /// Set the partition of the chemical system using a formatted string.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    auto setPartition(std::string partition) -> void;

    /// Solve the kinetic path problem.
    /// @param state The initial state of the chemical system
    /// @param t0 The initial time of the kinetic path
    /// @param t1 The final time of the kinetic path
    /// @param units The time units of `t0` and `t1` (e.g., `s`, `minute`, `day`, `year`, etc.).
    auto solve(ChemicalState& state, double t0, double t1, std::string units = "s") -> void;

    /// Return a ChemicalPlot instance.
    /// The returned ChemicalOutput instance must be properly configured
    /// before the method EquilibriumPath::solve is called.
    /// Changes in this ChemicalOutput instance are observed by the
    /// EquilibriumPath object.
    auto output() -> ChemicalOutput;

    /// Return a ChemicalPlot instance.
    /// The returned ChemicalPlot instance must be properly configured
    /// before the method EquilibriumPath::solve is called.
    /// Changes in this ChemicalPlot instance are observed by the
    /// EquilibriumPath object.
    auto plot() -> ChemicalPlot;

    /// Return a collection of ChemicalPlot instances.
    /// The returned ChemicalPlot instances must be properly configured
    /// before the method EquilibriumPath::solve is called.
    /// Changes in theses ChemicalPlot instances are observed by the
    /// EquilibriumPath object.
    auto plots(unsigned num) -> std::vector<ChemicalPlot>;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
