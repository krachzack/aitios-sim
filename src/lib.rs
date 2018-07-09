#[cfg_attr(test, macro_use)]
extern crate aitios_geom as geom;
extern crate aitios_scene as scene;
extern crate aitios_sampling as sampling;
extern crate aitios_spatial as spatial;
extern crate aitios_surf as surf;
extern crate rand;
#[macro_use]
extern crate log;

mod sim;
mod surfel_data;
mod surfel_rule;
mod ton;
mod tracer;

pub use sim::Simulation;
pub use ton::{TonSource, TonSourceBuilder};
pub use surfel_data::SurfelData;
pub use surfel_rule::SurfelRule;

#[cfg(feature = "export_tracer")]
pub use tracer::*;
