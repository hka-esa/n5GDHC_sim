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
