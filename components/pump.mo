/*

This file is part of n5GDHC_sim.

Copyright (c) 2025 Manuel Kollmar, Adrian Bürger, Markus Bohlayer, Marco Braun, Moritz Diehl.
Developed at HS Karlsruhe and IMTEK, University of Freiburg.
All rights reserved.

The BSD 3-Clause License

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

within n5GDHC.components;

model pump
  connectors.fluid_p outlet annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}})));
  connectors.fluid_p inlet annotation(
    Placement(transformation(extent = {{50, -10}, {70, 10}})));
  Modelica.Blocks.Interfaces.RealInput mdot annotation(
    Placement(visible = true, transformation(origin = {2, -22}, extent = {{-58, 52}, {-18, 92}}, rotation = 0), iconTransformation(origin = {1, 61}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
equation
  assert(abs(outlet.mdot - inlet.mdot) <= 1e-5, "Pump mass balance violated.", level = AssertionLevel.error);
  outlet.mdot = mdot;
  outlet.T = inlet.T;
  outlet.p = inlet.p;
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -80}, {80, 80}}), graphics = {Ellipse(lineThickness = 0.5, extent = {{-60, 60}, {60, -60}}), Line(points = {{-56, 20}, {20, 56}}, thickness = 0.5), Line(points = {{-56, -20}, {20, -56}}, thickness = 0.5)}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -80}, {80, 80}})));
end pump;
