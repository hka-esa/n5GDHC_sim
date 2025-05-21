within n5GDHC.components;

model hn_split
  connectors.fluid_p inlet annotation(
    Placement(transformation(extent = {{-50, -10}, {-30, 10}}), iconTransformation(extent = {{-50, -10}, {-30, 10}})));
  connectors.fluid_p outlet_1 annotation(
    Placement(transformation(extent = {{-10, 30}, {10, 50}}), iconTransformation(extent = {{-10, 30}, {10, 50}})));
  connectors.fluid_p outlet_2 annotation(
    Placement(transformation(extent = {{30, -10}, {50, 10}}), iconTransformation(extent = {{30, -10}, {50, 10}})));
equation
  assert(outlet_1.mdot >= (-1e-5), "hn_tube_s O1 negative massflow", level = AssertionLevel.error);
  assert(outlet_2.mdot >= (-1e-5), "hn_tube_s O2 negative massflow", level = AssertionLevel.error);
  outlet_1.T = inlet.T;
  outlet_2.T = inlet.T;
  outlet_1.mdot + outlet_2.mdot = inlet.mdot;
  outlet_1.p = inlet.p;
  outlet_2.p = inlet.p;
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-60, -60}, {60, 60}}), graphics = {Rectangle( lineThickness = 0.5, extent = {{-40, 40}, {40, -40}}), Line(points = {{22, 0}, {16, 6}, {16, -6}, {22, 0}}, thickness = 0.5), Line(points = {{-22, 0}, {16, 0}}, thickness = 0.5), Line(origin = {0, 21}, rotation = 90,points = {{3, 0}, {-3, 6}, {-3, -6}, {3, 0}}, thickness = 0.5), Line(points = {{0, 18}, {0, 0}}, thickness = 0.5)}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-60, -60}, {60, 60}})));
end hn_split;
