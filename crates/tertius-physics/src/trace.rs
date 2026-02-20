//! Derivation trace model and renderers.

/// A single symbolic rewrite/explanation step.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DerivationStep {
    /// Human-readable rule identifier.
    pub rule_applied: String,
    /// Expression before rule application.
    pub before: String,
    /// Expression after rule application.
    pub after: String,
    /// Optional guard conditions.
    pub conditions: Vec<String>,
}

/// A complete derivation for a solved problem.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DerivationTrace {
    /// Title for textbook output.
    pub title: String,
    /// Ordered derivation steps.
    pub steps: Vec<DerivationStep>,
    /// Final result string.
    pub final_result: String,
}

impl DerivationTrace {
    /// Creates an empty trace.
    #[must_use]
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            steps: Vec::new(),
            final_result: String::new(),
        }
    }

    /// Appends a new derivation step.
    pub fn push_step(
        &mut self,
        rule_applied: impl Into<String>,
        before: impl Into<String>,
        after: impl Into<String>,
        conditions: Vec<String>,
    ) {
        self.steps.push(DerivationStep {
            rule_applied: rule_applied.into(),
            before: before.into(),
            after: after.into(),
            conditions,
        });
    }

    /// Renders the derivation as textbook-style plain text.
    #[must_use]
    pub fn render_textbook(&self) -> String {
        let mut out = String::new();
        out.push_str(&self.title);
        out.push('\n');
        out.push('\n');

        for (idx, step) in self.steps.iter().enumerate() {
            out.push_str(&format!("Step {} ({})\n", idx + 1, step.rule_applied));
            out.push_str(&format!("Before: {}\n", step.before));
            out.push_str(&format!("After: {}\n", step.after));
            if !step.conditions.is_empty() {
                out.push_str(&format!("Conditions: {}\n", step.conditions.join(", ")));
            }
            out.push('\n');
        }

        if !self.final_result.is_empty() {
            out.push_str(&format!("Final: {}", self.final_result));
        }
        out
    }

    /// Renders the derivation in LaTeX align format.
    #[must_use]
    pub fn render_latex(&self) -> String {
        let mut out = String::new();
        out.push_str("\\begin{align*}\n");
        out.push_str(&format!(
            "\\text{{{}}}\\\\\n",
            self.title.replace('{', "\\{").replace('}', "\\}")
        ));
        for (idx, step) in self.steps.iter().enumerate() {
            out.push_str(&format!(
                "\\text{{Step {}: }} {} &\\Rightarrow {}\\\\\n",
                idx + 1,
                step.before,
                step.after
            ));
            out.push_str(&format!(
                "\\text{{Rule: {}}}\\\\\n",
                step.rule_applied.replace('{', "\\{").replace('}', "\\}")
            ));
            if !step.conditions.is_empty() {
                out.push_str(&format!(
                    "\\text{{Conditions: {}}}\\\\\n",
                    step.conditions
                        .join(", ")
                        .replace('{', "\\{")
                        .replace('}', "\\}")
                ));
            }
        }
        if !self.final_result.is_empty() {
            out.push_str(&format!("\\text{{Final: }} {}\\\\\n", self.final_result));
        }
        out.push_str("\\end{align*}");
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_render_textbook_includes_rule_provenance() {
        let mut trace = DerivationTrace::new("First-Order ODE");
        trace.push_step(
            "Integrating factor",
            "y' + a y = b",
            "(mu y)' = b mu",
            vec!["a != 0".to_string()],
        );
        trace.final_result = "y(x)=b/a + C e^{-ax}".to_string();
        let rendered = trace.render_textbook();
        assert!(rendered.contains("Step 1"));
        assert!(rendered.contains("Integrating factor"));
        assert!(rendered.contains("a != 0"));
        assert!(rendered.contains("Final: y(x)=b/a + C e^{-ax}"));
    }

    #[test]
    fn test_render_latex_wraps_steps() {
        let mut trace = DerivationTrace::new("Heat Equation");
        trace.push_step(
            "Separation of variables",
            "u_t = alpha^2 u_xx",
            "X T' = alpha^2 X'' T",
            vec![],
        );
        trace.final_result = "u(x,t)=sum_n A_n sin(n pi x/L)e^{-...}".to_string();
        let latex = trace.render_latex();
        assert!(latex.contains("\\begin{align*}"));
        assert!(latex.contains("\\text{Rule: Separation of variables}"));
        assert!(latex.contains("\\text{Final: }"));
        assert!(latex.contains("\\end{align*}"));
    }
}
