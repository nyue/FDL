// CS 5643 p4 RigidFluid

#include "textures.inc"
#include "colors.inc"
#include "functions.inc"
#include "glass.inc"  

               
#declare r1 = seed(0);

#declare WaterMaterial = material {
 texture {
  pigment { color rgbf <0.98, 0.98, 0.98, 0.98> }
  normal{ waves .05 }
  finish{ ambient 0 reflection { 0.1, 1 fresnel } conserve_energy phong 1 phong_size 90 }}
 interior{ ior 1.33335 caustics 1.0 fade_power 1 dispersion 1.01 dispersion_samples 16 }
};

#declare Ground =
box { <-100.0, -2.0, -100.0>, <100.0, -1.0, 100.0>
   texture {
      pigment {
         color rgb <0.98, 0.83, 0.58>
      }
      scale 5
      translate 10
   }
}

sky_sphere {pigment {rgb <0.5,0.55,0.7>}}

object {Ground}

//----------------------- THE TABLE
#declare Pig_1 =
pigment {
   gradient z
   color_map {
      [0.00, rgb <0.01, 0.59, 0.81>]
      [0.70, rgb <0.01, 0.59, 0.81>]
      [0.70, rgb <0.98, 0.98, 0.87>]
      [1.00, rgb <0.98, 0.98, 0.87>]
   }
   frequency 4
}

#declare Pig_2 =
pigment {
   bozo
   color_map {
      [0.00, rgb <0.35, 0.58, 0.88>*1.0]
      [0.25, rgb <0.35, 0.58, 0.88>*1.1]
      [0.50, rgb <0.35, 0.58, 0.88>*0.9]
      [0.75, rgb <0.35, 0.58, 0.88>*1.0]
      [1.00, rgb <0.35, 0.58, 0.88>*0.8]
   }
   scale 0.1
}

pigment
{
   density_file df3 "density_export.df3"
   interpolate 1 
   // [PIGMENT_MODIFIERS...]
}

light_source {<-34,40,50> rgb <0.9, 0.9, 0.95>*0.8}
light_source {
  <-10,50,10> rgb <1.0, 1.0, 0.95>
  area_light
  <1,0,0>, <0,0,1>, 3, 3
  adaptive 1
  shadowless
}
