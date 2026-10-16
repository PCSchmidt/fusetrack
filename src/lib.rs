//! FuseTrack: multi-sensor fusion tracking for counter-UAS.
//!
//! v0.1 provides a scenario simulator and a single-target EKF fed by a
//! range-bearing radar. Multi-target association, IMM and additional sensors
//! are on the roadmap (see `docs/roadmap.md`).

pub mod ekf;
pub mod models;
pub mod sim;

pub use ekf::Ekf;
pub use models::{wrap_angle, Radar, RadarMeas, State};
