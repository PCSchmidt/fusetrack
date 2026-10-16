//! FuseTrack: multi-sensor fusion tracking for counter-UAS.
//!
//! Provides a multi-target scenario simulator (missed detections and
//! clutter), a single-target EKF, an IMM filter for manoeuvring targets, a GNN
//! tracker with gating and track management, and OSPA scoring. JPDA and
//! additional sensors are on the roadmap (see `docs/roadmap.md`).

pub mod assignment;
pub mod ekf;
pub mod eval;
pub mod imm;
pub mod metrics;
pub mod models;
pub mod sim;
pub mod tracker;

pub use ekf::Ekf;
pub use imm::{Imm, ImmConfig, MotionModel};
pub use models::{wrap_angle, Radar, RadarMeas, State};
pub use tracker::{GnnTracker, ScoreConfig, TrackerConfig};
