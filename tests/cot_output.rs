//! A confirmed track becomes a valid unknown-UAS CoT event for TAK.

#![cfg(feature = "cot")]

use fusetrack::cot::{track_to_cot, Origin};
use fusetrack::ekf::Ekf;
use fusetrack::tracker::{Track, TrackFilter, TrackStatus};
use nalgebra::{Matrix4, Vector4};

fn track_at(id: u64, status: TrackStatus, state: Vector4<f64>) -> Track {
    let ekf = Ekf::new(state, Matrix4::identity() * 100.0, 1.0);
    Track {
        id,
        filter: TrackFilter::Single(ekf),
        status,
        age: 10,
        hits: 8,
        misses: 0,
        score: 0.0,
    }
}

#[test]
fn confirmed_track_becomes_unknown_uas_cot() {
    // Origin near DCA; track 1000 m east, 2000 m north, heading north-east.
    let origin = Origin {
        lat: 38.8512,
        lon: -77.0402,
        hae: 50.0,
    };
    let track = track_at(
        7,
        TrackStatus::Confirmed,
        Vector4::new(1000.0, 2000.0, 15.0, 15.0),
    );

    let event = track_to_cot(&track, &origin, 1_700_000_000_000_000, 30)
        .expect("a confirmed track is published");
    let xml = event.to_xml();

    // Correctly typed as an unknown UAS air track for TAK.
    assert!(xml.contains("type=\"a-u-A-M-F-Q\""));
    assert!(xml.contains("uid=\"fusetrack-7\""));
    // 2000 m north of 38.8512 is about +0.018 deg latitude; east of the origin.
    assert!(xml.contains("lat=\"38.86"), "unexpected point: {xml}");
    assert!(xml.contains("lon=\"-77.02"), "unexpected point: {xml}");
    // Heading north-east is 45 deg true at ~21 m/s.
    assert!(xml.contains("course=\"45.0\""), "unexpected track: {xml}");
    assert!(xml.contains("speed=\"21.21\""), "unexpected track: {xml}");
}

#[test]
fn tentative_track_is_not_published() {
    let origin = Origin {
        lat: 38.8512,
        lon: -77.0402,
        hae: 50.0,
    };
    let track = track_at(8, TrackStatus::Tentative, Vector4::new(0.0, 0.0, 0.0, 0.0));
    assert!(track_to_cot(&track, &origin, 1_700_000_000_000_000, 30).is_none());
}
