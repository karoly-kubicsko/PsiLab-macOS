//
//  LaTeXText.swift
//  MacQC
//
//  ΨLab 9.0: names of custom molecules written in LaTeX (the comment line of an .xyz file or the
//  title of a Gaussian input), e.g. "CO$_2$", "NH$_4^+$", "SO$_4^{2-}$", "$\alpha$-D-glucose",
//  "\ce{H2SO4}". Supported: math mode $…$ (letters italic, as in LaTeX; "-" a minus sign),
//  subscripts and superscripts (_x, _{…}, ^x, ^{…}, also outside math mode), \mathrm, \text,
//  \textrm, \mathbf, \textbf, \mathit, \textit, \emph, \mathsf; Greek letters; common symbols
//  (\cdot, \times, \pm, \to, \rightarrow, \leftrightarrow, \rightleftharpoons, \circ, \degree,
//  \AA, \infty, \prime, \ldots); spacing (\, \; \quad ~); escaped characters (\$ \_ \% \& \{ \});
//  the mhchem command \ce{…} (digits after an element or a bracket become subscripts, a final
//  charge a superscript, -> an arrow). Unknown commands are shown as written.
//
//  LaTeX.runs parses the text into runs (base line, subscript, superscript; italic, bold);
//  LaTeX.text draws them as SwiftUI Text (smaller, raised or lowered), LaTeX.plain gives a
//  plain-text form with Unicode subscripts and superscripts (menus, window titles, exports).
//

import SwiftUI
import AppKit

nonisolated enum LaTeX {
    enum Position: Equatable, Sendable {
        case base, sub, sup
    }

    struct Run: Equatable, Sendable {
        var text: String
        var position: Position
        var italic: Bool
        var bold: Bool
    }

    /// Font style of a group.
    private struct Style {
        var position: Position = .base
        var math = false        // math mode: letters italic, spaces ignored, "-" is a minus sign
        var upright = false     // \mathrm, \text: letters upright in math mode
        var italic = false      // \textit, \mathit, \emph
        var bold = false
        var text = false        // \text{…} inside math mode: spaces kept
    }

    private static let greek: [String: String] = [
        "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ϵ", "varepsilon": "ε", "zeta": "ζ",
        "eta": "η", "theta": "θ", "vartheta": "ϑ", "iota": "ι", "kappa": "κ", "lambda": "λ", "mu": "μ", "nu": "ν",
        "xi": "ξ", "pi": "π", "varpi": "ϖ", "rho": "ρ", "varrho": "ϱ", "sigma": "σ", "varsigma": "ς", "tau": "τ",
        "upsilon": "υ", "phi": "ϕ", "varphi": "φ", "chi": "χ", "psi": "ψ", "omega": "ω",
        "Gamma": "Γ", "Delta": "Δ", "Theta": "Θ", "Lambda": "Λ", "Xi": "Ξ", "Pi": "Π", "Sigma": "Σ",
        "Upsilon": "Υ", "Phi": "Φ", "Psi": "Ψ", "Omega": "Ω",
    ]

    private static let symbols: [String: String] = [
        "cdot": "·", "times": "×", "pm": "±", "mp": "∓", "to": "→", "rightarrow": "→", "leftarrow": "←",
        "leftrightarrow": "↔", "rightleftharpoons": "⇌", "Rightarrow": "⇒", "longrightarrow": "⟶",
        "uparrow": "↑", "downarrow": "↓", "circ": "∘", "degree": "°", "AA": "Å", "infty": "∞", "prime": "′",
        "ldots": "…", "dots": "…", "cdots": "⋯", "equiv": "≡", "approx": "≈", "neq": "≠", "ne": "≠", "leq": "≤",
        "le": "≤", "geq": "≥", "ge": "≥", "partial": "∂", "nabla": "∇", "hbar": "ℏ", "ell": "ℓ", "bullet": "•",
        "ddagger": "‡", "dagger": "†", "ast": "∗", "star": "⋆", "langle": "⟨", "rangle": "⟩", "S": "§",
        "textdegree": "°", "textendash": "–", "textemdash": "—", "minus": "−", "deg": "°",
    ]

    private static let fontCommands: Set<String> = [
        "mathrm", "text", "textrm", "mathbf", "textbf", "mathit", "textit", "emph", "mathsf", "textsf",
        "mathnormal", "operatorname", "textup", "mathup",
    ]

    // MARK: Parsing

    /// The runs of a LaTeX text, adjacent runs of the same style merged.
    static func runs(_ source: String) -> [Run] {
        var parser = Parser(chars: Array(source))
        parser.parseSequence(Style(), until: nil)
        return merged(parser.output)
    }

    private static func merged(_ runs: [Run]) -> [Run] {
        var out: [Run] = []
        for r in runs where !r.text.isEmpty {
            if let last = out.last, last.position == r.position, last.italic == r.italic, last.bold == r.bold {
                out[out.count - 1].text += r.text
            } else {
                out.append(r)
            }
        }
        return out
    }

    private struct Parser {
        let chars: [Character]
        var i = 0
        var output: [Run] = []

        init(chars: [Character]) {
            self.chars = chars
        }

        mutating func emit(_ s: String, _ style: Style, letterLike: Bool = false) {
            guard !s.isEmpty else { return }
            let italic = style.italic || (style.math && !style.upright && letterLike)
            output.append(Run(text: s, position: style.position, italic: italic, bold: style.bold))
        }

        /// Parses until `end` ("}" of a group or "$" of math mode) or the end of the text.
        mutating func parseSequence(_ style0: Style, until end: Character?) {
            var style = style0
            while i < chars.count {
                let c = chars[i]
                if let e = end, c == e {
                    i += 1
                    return
                }
                switch c {
                case "$":
                    i += 1
                    if i < chars.count && chars[i] == "$" { i += 1 }  // $$ … $$ like $ … $
                    if style.math {
                        // closing $ without an open group: leave math mode
                        style.math = false
                    } else {
                        var m = style
                        m.math = true
                        m.upright = false
                        m.text = false
                        parseSequence(m, until: "$")
                    }
                case "{":
                    i += 1
                    parseSequence(style, until: "}")
                case "}":
                    i += 1  // unbalanced: ignored
                case "_", "^":
                    i += 1
                    var s = style
                    s.position = c == "_" ? .sub : .sup
                    parseArgument(s)
                case "\\":
                    parseCommand(&style)
                case "~":
                    i += 1
                    emit("\u{00A0}", style)
                case " ", "\t", "\n":
                    i += 1
                    if !style.math || style.text { emit(" ", style) }
                case "-":
                    i += 1
                    if style.math && !style.text {
                        emit("−", style)
                    } else if i + 1 < chars.count && chars[i] == "-" && chars[i + 1] == "-" {
                        i += 2
                        emit("—", style)
                    } else if i < chars.count && chars[i] == "-" {
                        i += 1
                        emit("–", style)
                    } else {
                        emit("-", style)
                    }
                case "'":
                    i += 1
                    emit(style.math && !style.text ? "′" : "'", style)
                default:
                    i += 1
                    emit(String(c), style, letterLike: c.isLetter)
                }
            }
        }

        /// The argument of _, ^ or a command: a {group}, a command or one character.
        mutating func parseArgument(_ style: Style) {
            while i < chars.count && (chars[i] == " " || chars[i] == "\t") { i += 1 }
            guard i < chars.count else { return }
            let c = chars[i]
            if c == "{" {
                i += 1
                parseSequence(style, until: "}")
            } else if c == "\\" {
                var s = style
                parseCommand(&s)
            } else if c == "-" && style.math && !style.text {
                i += 1
                emit("−", style)
            } else {
                i += 1
                emit(String(c), style, letterLike: c.isLetter)
            }
        }

        mutating func parseCommand(_ style: inout Style) {
            i += 1  // the backslash
            guard i < chars.count else {
                emit("\\", style)
                return
            }
            let c = chars[i]
            if !c.isLetter {
                i += 1
                switch c {
                case ",": emit("\u{2009}", style)
                case ";", ":", ">": emit("\u{2005}", style)
                case "!": break
                case " ": emit(" ", style)
                case "\\": emit(" ", style)
                default: emit(String(c), style)  // \$ \_ \% \& \{ \} \#
                }
                return
            }
            var name = ""
            while i < chars.count && chars[i].isLetter {
                name.append(chars[i])
                i += 1
            }
            if LaTeX.fontCommands.contains(name) {
                var s = style
                switch name {
                case "mathrm", "textrm", "mathsf", "textsf", "operatorname", "textup", "mathup":
                    s.upright = true
                    s.italic = false
                case "text":
                    s.upright = true
                    s.italic = false
                    s.text = true
                case "mathbf", "textbf":
                    s.bold = true
                    s.upright = true
                case "mathit", "textit", "emph":
                    s.italic = true
                default:
                    break
                }
                parseArgument(s)
                return
            }
            switch name {
            case "rm", "upshape":
                style.upright = true
                style.italic = false
                skipSpace()
            case "it", "itshape":
                style.italic = true
                skipSpace()
            case "bf", "bfseries":
                style.bold = true
                skipSpace()
            case "ce":
                var s = style
                s.upright = true
                s.italic = false
                s.text = true
                parseChemistry(s)
            case "quad":
                emit("\u{2003}", style)
            case "qquad":
                emit("\u{2003}\u{2003}", style)
            case "left", "right", "big", "Big", "bigl", "bigr", "displaystyle", "textstyle", "mathchoice":
                skipSpace()
            default:
                if let g = LaTeX.greek[name] {
                    skipSpace()
                    emit(g, style, letterLike: name.first?.isLowercase == true)
                } else if name == "circ" && style.position == .sup {
                    skipSpace()
                    emit("°", style)
                } else if let s = LaTeX.symbols[name] {
                    skipSpace()
                    emit(s, style)
                } else {
                    emit("\\" + name, style)
                }
            }
        }

        mutating func skipSpace() {
            while i < chars.count && chars[i] == " " { i += 1 }
        }

        /// mhchem \ce{…}: "H2SO4", "SO4^2-", "NH4+", "Fe(CN)6^3-", "A -> B", "2 H2O".
        mutating func parseChemistry(_ style: Style) {
            while i < chars.count && chars[i] == " " { i += 1 }
            guard i < chars.count && chars[i] == "{" else { return }
            i += 1
            var body: [Character] = []
            var depth = 0
            while i < chars.count {
                let c = chars[i]
                i += 1
                if c == "{" { depth += 1 }
                if c == "}" {
                    if depth == 0 { break }
                    depth -= 1
                }
                body.append(c)
            }
            var k = 0
            var previousIsAtom = false  // a letter, ")" or "]" before: digits are subscripts
            while k < body.count {
                let c = body[k]
                if c == "-" && k + 1 < body.count && body[k + 1] == ">" {
                    emit(" → ", style)
                    k += 2
                    previousIsAtom = false
                } else if c == "<" && k + 2 < body.count && body[k + 1] == "=" && body[k + 2] == ">" {
                    emit(" ⇌ ", style)
                    k += 3
                    previousIsAtom = false
                } else if c == "^" {
                    var s = style
                    s.position = .sup
                    k += 1
                    if k < body.count && body[k] == "{" {
                        k += 1
                        var t = ""
                        while k < body.count && body[k] != "}" {
                            t.append(body[k] == "-" ? "−" : body[k])
                            k += 1
                        }
                        k += 1
                        emit(t, s)
                    } else {
                        var t = ""
                        while k < body.count && (body[k].isNumber || body[k] == "+" || body[k] == "-") {
                            t.append(body[k] == "-" ? "−" : body[k])
                            k += 1
                        }
                        emit(t, s)
                    }
                    previousIsAtom = false
                } else if c == "_" {
                    var s = style
                    s.position = .sub
                    k += 1
                    var t = ""
                    if k < body.count && body[k] == "{" {
                        k += 1
                        while k < body.count && body[k] != "}" {
                            t.append(body[k])
                            k += 1
                        }
                        k += 1
                    } else if k < body.count {
                        t.append(body[k])
                        k += 1
                    }
                    emit(t, s)
                } else if c.isNumber && previousIsAtom {
                    var t = ""
                    while k < body.count && body[k].isNumber {
                        t.append(body[k])
                        k += 1
                    }
                    var s = style
                    s.position = .sub
                    emit(t, s)
                } else if (c == "+" || c == "-") && previousIsAtom
                            && (k + 1 == body.count || body[k + 1] == " ") {
                    var s = style
                    s.position = .sup
                    emit(c == "-" ? "−" : "+", s)
                    k += 1
                } else {
                    emit(String(c), style)
                    previousIsAtom = c.isLetter || c == ")" || c == "]"
                    k += 1
                }
            }
        }
    }

    // MARK: Output

    private static let subscripts: [Character: Character] = [
        "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉",
        "+": "₊", "-": "₋", "−": "₋", "=": "₌", "(": "₍", ")": "₎", "a": "ₐ", "e": "ₑ", "o": "ₒ", "x": "ₓ",
        "h": "ₕ", "k": "ₖ", "l": "ₗ", "m": "ₘ", "n": "ₙ", "p": "ₚ", "s": "ₛ", "t": "ₜ", "i": "ᵢ", "r": "ᵣ",
        "u": "ᵤ", "v": "ᵥ", "j": "ⱼ", "β": "ᵦ", "γ": "ᵧ", "ρ": "ᵨ", "φ": "ᵩ", "χ": "ᵪ",
    ]

    private static let superscripts: [Character: Character] = [
        "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
        "+": "⁺", "-": "⁻", "−": "⁻", "=": "⁼", "(": "⁽", ")": "⁾", "n": "ⁿ", "i": "ⁱ", "a": "ᵃ", "b": "ᵇ",
        "c": "ᶜ", "d": "ᵈ", "e": "ᵉ", "f": "ᶠ", "g": "ᵍ", "h": "ʰ", "j": "ʲ", "k": "ᵏ", "l": "ˡ", "m": "ᵐ",
        "o": "ᵒ", "p": "ᵖ", "r": "ʳ", "s": "ˢ", "t": "ᵗ", "u": "ᵘ", "v": "ᵛ", "w": "ʷ", "x": "ˣ", "y": "ʸ",
        "z": "ᶻ", "A": "ᴬ", "B": "ᴮ", "D": "ᴰ", "E": "ᴱ", "G": "ᴳ", "H": "ᴴ", "I": "ᴵ", "J": "ᴶ", "K": "ᴷ",
        "L": "ᴸ", "M": "ᴹ", "N": "ᴺ", "O": "ᴼ", "P": "ᴾ", "R": "ᴿ", "T": "ᵀ", "U": "ᵁ", "V": "ⱽ", "W": "ᵂ",
        "∘": "°", "°": "°", "′": "′", "*": "*", "‡": "‡",
    ]

    /// Plain text with Unicode subscripts and superscripts ("CO$_2$" -> "CO₂", "NH$_4^+$" -> "NH₄⁺");
    /// characters without a Unicode form are written with _ or ^ ("x_{ab}" -> "x_ab").
    static func plain(_ source: String) -> String {
        var out = ""
        for r in runs(source) {
            switch r.position {
            case .base:
                out += r.text
            case .sub, .sup:
                let table = r.position == .sub ? subscripts : superscripts
                let mapped = r.text.map { table[$0] }
                if mapped.allSatisfy({ $0 != nil }) {
                    out += String(mapped.compactMap { $0 })
                } else {
                    out += (r.position == .sub ? "_" : "^") + r.text
                }
            }
        }
        return out
    }

    /// The runs as one SwiftUI Text: subscripts and superscripts at 70 % of `size`, lowered by
    /// 0.2 and raised by 0.38 of it.
    static func text(_ source: String, size: CGFloat = 13, weight: Font.Weight = .regular) -> Text {
        var result = Text(verbatim: "")
        for r in runs(source) {
            let small = r.position != .base
            var t = Text(verbatim: r.text)
                .font(.system(size: small ? size * 0.7 : size, weight: r.bold ? .bold : weight))
            if r.italic { t = t.italic() }
            switch r.position {
            case .base: break
            case .sub: t = t.baselineOffset(-size * 0.2)
            case .sup: t = t.baselineOffset(size * 0.38)
            }
            result = result + t
        }
        return result
    }
}

/// A LaTeX text drawn with subscripts and superscripts (molecule names of custom molecules).
struct LaTeXText: View {
    let source: String
    var size: CGFloat = 13
    var weight: Font.Weight = .regular

    init(_ source: String, size: CGFloat = 13, weight: Font.Weight = .regular) {
        self.source = source
        self.size = size
        self.weight = weight
    }

    var body: some View {
        LaTeX.text(source, size: size, weight: weight)
            .accessibilityLabel(LaTeX.plain(source))
    }
}

/// 20.0: the note under the molecule picker (electrons, term symbol, point group, geometry), e.g.
/// "1 electron, X ^2Σ_g^+, R = 1.0569 Å (2.00 bohr)" or "…, C_{3v}, …": caption size, secondary
/// colour, with the subscripts and superscripts of LaTeX (Σ_g^+, C_{3v} as one subscript).
struct MoleculeCaption: View {
    let source: String

    init(_ source: String) {
        self.source = source
    }

    var body: some View {
        LaTeX.text(source, size: NSFont.preferredFont(forTextStyle: .caption1).pointSize)
            .foregroundStyle(.secondary)
            .fixedSize(horizontal: false, vertical: true)
            .accessibilityLabel(LaTeX.plain(source))
    }
}
