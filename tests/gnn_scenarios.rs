//! End-to-end regression tests for the GNN tracker on simulated scenarios.
//!
//! Bounds sit comfortably above the values reported in the README so that
//! they catch real regressions without being brittle.

use fusetrack::eval::run_gnn_mean;
use fusetrack::sim::ScenarioConfig;
use fusetrack::{Radar, TrackerConfig};
use nalgebra::Vector2;

fn scenario(pd: f64, clutter_rate: f64) -> ScenarioConfig {
    ScenarioConfig {
        radar: Radar {
            pos: Vector2::new(0.0, 0.0),
            sigma_r: 15.0,
            sigma_b: 0.003,
        },
        dt: 1.0,
        steps: 100,
        n_targets: 5,
        sigma_accel: 1.0,
        pd,
        clutter_rate,
        r_max: 15000.0,
    }
}

#[test]
fn clean_scenario_tracks_every_target_with_small_error() {
    let s = run_gnn_mean(&scenario(1.0, 0.0), &TrackerConfig::default(), 5);
    println!("{s:?}");
    assert!(s.mean_ospa < 20.0, "OSPA {}", s.mean_ospa);
    assert!((s.mean_tracks - s.mean_targets).abs() < 0.1);
}

#[test]
fn light_clutter_and_missed_detections_are_handled() {
    let s = run_gnn_mean(&scenario(0.9, 10.0), &TrackerConfig::default(), 5);
    println!("{s:?}");
    assert!(s.mean_ospa < 25.0, "OSPA {}", s.mean_ospa);
    assert!((s.mean_tracks - s.mean_targets).abs() < 0.5);
}

#[test]
fn more_clutter_is_worse() {
    let light = run_gnn_mean(&scenario(0.9, 10.0), &TrackerConfig::default(), 5);
    let heavy = run_gnn_mean(&scenario(0.9, 50.0), &TrackerConfig::default(), 5);
    assert!(heavy.mean_ospa > light.mean_ospa);
}
