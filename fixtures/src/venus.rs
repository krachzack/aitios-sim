use asset::obj;
use sim::{Simulation, TonSourceBuilder, SurfelRule, SurfelData, Tracer};
use surf::{SurfaceBuilder, SurfelSampling};
use scene::{Entity, Mesh};

pub fn entitites() -> Vec<Entity> {
    obj::load("test-scenes/venus.obj")
            .unwrap()
}

pub fn make_tracer() -> Tracer {
    let ents = entitites();
    let tris = ents.iter().flat_map(|e| e.mesh.triangles());
    Tracer::new(tris)
}

/// Creates a moderately complex simulation with a scene  from a ~1.7MB OBJ (30832 triangles)
/// of a simplified version of the venus the milo standing on a flattened out icosphere.
///
/// Venus the Milo was obtained from Scan the world.
/// https://de.wikipedia.org/wiki/Datei:Scan_the_World_-_Venus_de_Milo.stl
pub fn make_simulation(emission_count: usize) -> Simulation {
    let scene = entitites();

    let source = &obj::load("test-scenes/buddha-scene-ton-source-mesh/buddha-scene-ton-source-sun.obj")
        .unwrap()[0];

    let surfel_prototype = SurfelData {
        entity_idx: 0,
        delta_straight: 0.1,
        delta_parabolic: 0.1,
        delta_flow: 0.1,
        substances: vec![0.0, 0.0],
        deposition_rates: vec![1.0, 1.0],
        rules: vec![]
    };

    let scene_triangles = scene.iter().flat_map(|e| e.mesh.triangles());
    let surface = SurfaceBuilder::new()
        .sampling(SurfelSampling::MinimumDistance(0.01))
        .sample_triangles(scene_triangles, &surfel_prototype)
        .build();

    let sources = vec![
        TonSourceBuilder::new()
            .emission_count(emission_count)
            .entity_shaped(&source, true)
            .p_straight(0.3)
            .p_parabolic(0.3)
            .p_flow(0.0)
            .substances(&vec![1.0, 0.0])
            .pickup_rates(vec![1.0, 1.0])
            .parabola_height(0.1)
            .flow_distance(0.05)
            .build()
    ];

    let scene_triangles = scene.iter().flat_map(|e| e.mesh.triangles());
    let rules = vec![SurfelRule::Transfer {
        source_substance_idx: 0,
        target_substance_idx: 1,
        factor: 0.1
    }];

    Simulation::new(sources, scene_triangles, surface, rules)
}
