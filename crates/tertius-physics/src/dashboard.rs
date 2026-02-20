//! Public physics dashboard reporting helpers.

use std::fs;
use std::path::Path;

/// Correctness metrics for a physics corpus run.
#[derive(Clone, Debug)]
pub struct CorrectnessReport {
    /// Number of solved items.
    pub solved: usize,
    /// Total corpus size.
    pub total: usize,
    /// Number of regression failures.
    pub regressions: usize,
}

/// Performance metrics for physics workflows.
#[derive(Clone, Debug)]
pub struct PerformanceReport {
    /// Median runtime in milliseconds.
    pub median_ms: f64,
    /// P95 runtime in milliseconds.
    pub p95_ms: f64,
}

/// Combined dashboard snapshot.
#[derive(Clone, Debug)]
pub struct PhysicsDashboard {
    /// RFC3339-like timestamp.
    pub generated_at: String,
    /// Correctness section.
    pub correctness: CorrectnessReport,
    /// Performance section.
    pub performance: PerformanceReport,
}

impl PhysicsDashboard {
    /// Correctness score in `[0,1]`.
    #[must_use]
    pub fn correctness_score(&self) -> f64 {
        if self.correctness.total == 0 {
            return 0.0;
        }
        self.correctness.solved as f64 / self.correctness.total as f64
    }

    /// Renders a public Markdown dashboard.
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let score_pct = self.correctness_score() * 100.0;
        format!(
            "# Tertius Physics Dashboard\n\n\
Generated: `{}`\n\n\
## Correctness\n\
- Score: {:.2}%\n\
- Solved: {}/{}\n\
- Regressions: {}\n\n\
## Performance\n\
- Median: {:.2} ms\n\
- P95: {:.2} ms\n",
            self.generated_at,
            score_pct,
            self.correctness.solved,
            self.correctness.total,
            self.correctness.regressions,
            self.performance.median_ms,
            self.performance.p95_ms
        )
    }

    /// Renders a JSON snapshot.
    #[must_use]
    pub fn render_json(&self) -> String {
        format!(
            "{{\"generated_at\":\"{}\",\"correctness\":{{\"solved\":{},\"total\":{},\"regressions\":{},\"score\":{:.6}}},\"performance\":{{\"median_ms\":{:.6},\"p95_ms\":{:.6}}}}}",
            escape_json(&self.generated_at),
            self.correctness.solved,
            self.correctness.total,
            self.correctness.regressions,
            self.correctness_score(),
            self.performance.median_ms,
            self.performance.p95_ms
        )
    }

    /// Writes Markdown dashboard to disk.
    pub fn write_markdown(&self, path: impl AsRef<Path>) -> std::io::Result<()> {
        fs::write(path, self.render_markdown())
    }

    /// Writes JSON snapshot to disk.
    pub fn write_json(&self, path: impl AsRef<Path>) -> std::io::Result<()> {
        fs::write(path, self.render_json())
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn tmp_file(name: &str) -> String {
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        format!("/tmp/{name}_{ts}.txt")
    }

    #[test]
    fn test_render_markdown_contains_main_metrics() {
        let dashboard = PhysicsDashboard {
            generated_at: "2026-02-20T00:00:00Z".to_string(),
            correctness: CorrectnessReport {
                solved: 180,
                total: 200,
                regressions: 3,
            },
            performance: PerformanceReport {
                median_ms: 120.0,
                p95_ms: 240.0,
            },
        };
        let md = dashboard.render_markdown();
        assert!(md.contains("90.00%"));
        assert!(md.contains("180/200"));
        assert!(md.contains("Median"));
        assert!(md.contains("120"));
    }

    #[test]
    fn test_render_json_contains_keys() {
        let dashboard = PhysicsDashboard {
            generated_at: "2026-02-20T00:00:00Z".to_string(),
            correctness: CorrectnessReport {
                solved: 9,
                total: 10,
                regressions: 1,
            },
            performance: PerformanceReport {
                median_ms: 100.0,
                p95_ms: 190.0,
            },
        };
        let json = dashboard.render_json();
        assert!(json.contains("\"generated_at\""));
        assert!(json.contains("\"correctness\""));
        assert!(json.contains("\"performance\""));
    }

    #[test]
    fn test_writers_emit_files() {
        let dashboard = PhysicsDashboard {
            generated_at: "2026-02-20T00:00:00Z".to_string(),
            correctness: CorrectnessReport {
                solved: 1,
                total: 1,
                regressions: 0,
            },
            performance: PerformanceReport {
                median_ms: 10.0,
                p95_ms: 20.0,
            },
        };
        let md_path = tmp_file("physics_dashboard");
        let json_path = tmp_file("physics_dashboard_json");
        dashboard.write_markdown(&md_path).unwrap();
        dashboard.write_json(&json_path).unwrap();
        let md = fs::read_to_string(&md_path).unwrap();
        let json = fs::read_to_string(&json_path).unwrap();
        assert!(!md.is_empty());
        assert!(!json.is_empty());
    }
}
