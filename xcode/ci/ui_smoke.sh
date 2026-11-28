#!/bin/bash
# UI smoke test for the built ΨLab.app on a GitHub macOS runner.
# Launches the app, screenshots each tab, dumps the visible text of every tab to the log,
# runs the potential-curve scan and checks the Edit/File menu commands; then relaunches with
# launch arguments for H2O (RHF), LiH (Full CI), an H2 potential curve, the refused Full CI
# of H2O, a long H2O calculation (warning) and the coupled-cluster methods of ΨLab 4.0
# (open method menu, H2O CCSD(T), LiH CCSDT, refused H2O CCSDTQ), and the job types of ΨLab 5.0
# (open job menu, H2O RHF optimisation, LiH CCSD(T) frequencies and thermochemistry, H2O
# opt+freq+SP with a CCSD(T)/aug-cc-pVTZ single point, the Output tab while a job runs, Cancel,
# the refused CCSDTQ frequencies of H2O), ΨLab 6.0 (the molecule of each optimisation step,
# played step by step; the animated normal modes of LiH and H2O: two frames 0.3 s apart must
# differ, and a few screenshots are printed into the log as small JPEGs), and ΨLab 7.0 (HeH+ with
# Full CI/aug-cc-pVTZ opt+freq, the relaxed torsion scan of H2O2 with its points played, the
# inversion transition state of NH3 and the trans transition state of H2O2 with CalcAll, the
# scan settings in the sidebar, the refused TS search of LiH) and ΨLab 8.0 (the orbitals:
# isosurfaces, a contour plot, Full CI natural orbitals, the Szabo-Ostlund orbitals), ΨLab 9.0 (custom
# molecules from XYZ and Gaussian input text, the kept optimised geometry, the molecule menu, the
# file panel, refusals), ΨLab 10.0 (the analytic CCSD Hessian of H2O, the Hessian setting), ΨLab 12.0
# (the 3D molecule viewer: sticks with atom numbers, spheres on black) and ΨLab 13.0 (an output file of
# psilab opened again: every tab and the orbitals from its archive). The screen of the
# runner is 1024 pt wide: the screenshots show that the six tabs fit at the minimum window width.
# Never fails the job: it only gathers evidence (screenshots + text) for review.
#
# Usage: ui_smoke.sh path/to/ΨLab.app output-dir

APP="$1"
OUT="$2"
mkdir -p "$OUT"

# Reads the text of the front window; optionally clicks a tab/button first.
cat > "$OUT/tab.applescript" <<'APPLESCRIPT'
on run argv
    set targetName to item 1 of argv
    set logText to ""
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            if targetName is "__LAYOUT__" then
                -- window frame and the frames of the tab buttons, segmented pickers and menus
                set w to window 1
                set wp to position of w
                set ws to size of w
                set logText to "WINDOW x " & (item 1 of wp) & " y " & (item 2 of wp) & " width " & (item 1 of ws) & " height " & (item 2 of ws) & linefeed
                set els to entire contents of w
                repeat with e in els
                    try
                        set r to role of e
                        if r is "AXRadioButton" or r is "AXTabGroup" or r is "AXMenuButton" or r is "AXPopUpButton" then
                            set n to ""
                            set d to ""
                            try
                                set n to name of e
                            end try
                            try
                                set d to description of e
                            end try
                            set p to position of e
                            set z to size of e
                            set logText to logText & r & " | " & n & " | " & d & " | x " & (item 1 of p) & " y " & (item 2 of p) & " w " & (item 1 of z) & " h " & (item 2 of z) & linefeed
                        end if
                    end try
                end repeat
                return logText
            end if
            if targetName is "__TREE__" then
                -- accessibility tree (everything except static text), to find how tabs are exposed
                set els to entire contents of window 1
                repeat with e in els
                    try
                        set r to role of e
                        if r is not "AXStaticText" then
                            set n to ""
                            set d to ""
                            try
                                set n to name of e
                            end try
                            try
                                set d to description of e
                            end try
                            set logText to logText & r & " | " & n & " | " & d & linefeed
                        end if
                    end try
                end repeat
                return logText
            end if
            if targetName is "__SCAN__" then
                -- the Run scan button has no accessible name: click the first AXButton
                -- that follows the "Step (Å)" pop-up in the tree
                set els to entire contents of window 1
                set seenStep to false
                repeat with e in els
                    try
                        set r to role of e
                        set n to ""
                        try
                            set n to name of e
                        end try
                        if r is "AXPopUpButton" and n is "Step (Å)" then
                            set seenStep to true
                        else if seenStep and r is "AXButton" then
                            click e
                            return "CLICKED button after Step pop-up"
                        end if
                    end try
                end repeat
                return "NOT FOUND: button after Step pop-up"
            end if
            if targetName is not "" then
                set clicked to false
                set els to entire contents of window 1
                repeat with e in els
                    if clicked then exit repeat
                    try
                        set r to role of e
                        if r is not "AXStaticText" and r is not "AXImage" and r is not "AXGroup" and r is not "AXScrollArea" then
                            set n to ""
                            set d to ""
                            set t to ""
                            try
                                set n to name of e
                            end try
                            try
                                set d to description of e
                            end try
                            try
                                set t to title of e
                            end try
                            if n is targetName or d is targetName or t is targetName then
                                click e
                                set clicked to true
                                set logText to logText & "CLICKED " & r & " | " & n & " | " & d & linefeed
                            end if
                        end if
                    end try
                end repeat
                if not clicked then set logText to logText & "NOT FOUND: " & targetName & linefeed
            end if
            delay 2
            set els to entire contents of window 1
            repeat with e in els
                try
                    if (role of e) is "AXStaticText" then
                        set logText to logText & (value of e as text) & linefeed
                    end if
                end try
            end repeat
            return logText
        end tell
    end tell
end run
APPLESCRIPT

# macOS has no `timeout`: run a command and kill it after N seconds.
with_timeout() {
    local secs="$1"; shift
    "$@" &
    local pid=$!
    ( sleep "$secs"; kill "$pid" 2>/dev/null ) > /dev/null 2>&1 &
    local watchdog=$!
    wait "$pid"
    local rc=$?
    kill "$watchdog" 2>/dev/null
    return $rc
}

# Screenshots go into the ui-smoke artifact only (PNG and a smaller JPEG); the log keeps the
# text of every view and the layout numbers, so it stays small enough to read through the API.
shot() {
    screencapture -x "$OUT/$1.png" || { echo "screencapture failed for $1"; return; }
    sips -s format jpeg -s formatOptions 60 -Z 1100 "$OUT/$1.png" --out "$OUT/$1.jpg" > /dev/null 2>&1 || true
    echo "screenshot: $1.png"
}

# ΨLab 6.0: two captures of the window $2 seconds apart (default 0.3); different files mean
# that something moved.
frames() {   # $1 = label, $2 = seconds between the captures
    screencapture -x -R0,60,1024,650 "$OUT/$1-frame1.png" || return 0
    sleep "${2:-0.3}"
    screencapture -x -R0,60,1024,650 "$OUT/$1-frame2.png" || return 0
    if cmp -s "$OUT/$1-frame1.png" "$OUT/$1-frame2.png"; then
        echo "FRAMES $1: identical (nothing moved)"
    else
        echo "FRAMES $1: different (animated)"
    fi
}

# A screenshot as a small JPEG, base64-encoded into the log for review.
print_shot() {   # $1 = screenshot name
    sips -s format jpeg -s formatOptions 45 -Z 760 "$OUT/$1.png" --out "$OUT/$1-log.jpg" > /dev/null 2>&1 || return 0
    echo "=== BEGIN $1.jpg ==="
    base64 -i "$OUT/$1-log.jpg" | fold -w 4000
    echo "=== END ==="
}

# Screen size and the layout of the window (do the six tabs fit next to the sidebar?).
layout() {   # $1 = label
    echo "::group::layout: $1"
    osascript -e 'tell application "Finder" to get bounds of window of desktop' 2>/dev/null | sed 's/^/SCREEN bounds /'
    with_timeout 60 osascript "$OUT/tab.applescript" "__LAYOUT__" | tee "$OUT/layout-$1.txt"
    echo "::endgroup::"
}

dump_tab() {   # $1 = file prefix, $2 = element to click ("" = none)
    echo "::group::UI text: $1"
    with_timeout 90 osascript "$OUT/tab.applescript" "$2" | tee "$OUT/$1.txt"
    echo "::endgroup::"
    shot "$1"
}

# Selects a tab or page through View > <item> (unambiguous: the sidebar and the Methods tab
# also have "Full CI" segments), then dumps the visible text and takes a screenshot.
cat > "$OUT/view-menu.applescript" <<'APPLESCRIPT'
on run argv
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            click menu item (item 1 of argv) of menu 1 of menu bar item "View" of menu bar 1
        end tell
    end tell
    return "VIEW MENU " & (item 1 of argv)
end run
APPLESCRIPT
# 21.0: scrolls the front window so that the static text $1 is visible (for the charts below the
# first screen of a tab): AXScrollToVisible on the text, else every vertical scroll bar to the end.
cat > "$OUT/scroll-to.applescript" <<'APPLESCRIPT'
on run argv
    set target to item 1 of argv
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            set els to entire contents of window 1
            repeat with e in els
                try
                    if (role of e) is "AXStaticText" and (value of e as text) is target then
                        perform action "AXScrollToVisible" of e
                        return "SCROLLED (AXScrollToVisible) to " & target
                    end if
                end try
            end repeat
            set moved to 0
            repeat with e in els
                try
                    if (role of e) is "AXScrollArea" then
                        repeat with b in (scroll bars of e)
                            try
                                if (orientation of b) is "AXVerticalOrientation" then
                                    set value of b to 1.0
                                    set moved to moved + 1
                                end if
                            end try
                        end repeat
                    end if
                end try
            end repeat
            return "SCROLL BARS set to the end: " & moved & " (" & target & " not scrollable to)"
        end tell
    end tell
end run
APPLESCRIPT

scroll_shot() {   # $1 = screenshot name, $2 = static text to scroll to
    with_timeout 90 osascript "$OUT/scroll-to.applescript" "$2" || echo "could not scroll to $2"
    sleep 2
    shot "$1"
    print_shot "$1"
}

view_tab() {   # $1 = file prefix, $2 = View menu item
    with_timeout 20 osascript "$OUT/view-menu.applescript" "$2" || echo "could not select $2"
    sleep 2
    dump_tab "$1" ""
}

# 17.0: without launch arguments ΨLab opens with the H2 Full CI example; HeH+ in STO-nG is the
# Szabo-Ostlund calculation only with -psilabSzaboOstlund YES (off by default)
echo "== launching $APP (HeH+ STO-3G, Szabo-Ostlund)"
open -n "$APP" --args -psilabMolecule HeH+ -psilabBasis STO-3G -psilabSzaboOstlund YES
sleep 15
pgrep -x "ΨLab" >/dev/null && echo "ΨLab process is running" || echo "ΨLab process NOT running"

shot "00-launch"

# 20.0: the sections of the newest version first, so that the time limit of the job cannot cut them
# off; then the app is launched again as above for the sections of the earlier versions
# --- ΨLab 20.0: H3+ run after the start-up example (the sidebar's menus, Cmd-R): the 3D viewers of
# Summary, Optimisation and the normal modes must show the triangle (19.0: only the background) ---
cat > "$OUT/popup.applescript" <<'APPLESCRIPT'
on run argv
    set popName to item 1 of argv
    set itemName to item 2 of argv
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            set els to entire contents of window 1
            repeat with e in els
                try
                    if (role of e) is "AXPopUpButton" then
                        set n to ""
                        set d to ""
                        try
                            set n to name of e
                        end try
                        try
                            set d to description of e
                        end try
                        if n is popName or d is popName then
                            click e
                            delay 1
                            click menu item itemName of menu 1 of e
                            return "SELECTED " & itemName & " in " & popName
                        end if
                    end if
                end try
            end repeat
            return "NOT FOUND: " & popName
        end tell
    end tell
end run
APPLESCRIPT
echo "== 20.0: start-up example, then H3+ (Full CI/aug-cc-pVDZ opt+freq+sp) chosen in the sidebar"
pkill -x "ΨLab" || true
sleep 3
open -n "$APP" --args -appearanceChoice dark
sleep 20
view_tab "e0-example-optimisation" "Optimisation"
print_shot "e0-example-optimisation"
view_tab "e0-example-frequencies" "Frequencies"
view_tab "e0-example-summary" "Summary"
with_timeout 120 osascript "$OUT/popup.applescript" "Atom/Molecule" "H₃⁺" || echo "popup failed"
sleep 3
with_timeout 120 osascript "$OUT/popup.applescript" "Basis set" "aug-cc-pVDZ" || echo "popup failed"
sleep 3
with_timeout 120 osascript "$OUT/popup.applescript" "SP basis set" "aug-cc-pVDZ" || echo "popup failed"
sleep 3
shot "e1-h3-chosen"
print_shot "e1-h3-chosen"
cat > "$OUT/run.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        set frontmost to true
        keystroke "r" using command down
        delay 3
        try
            click button "Run" of sheet 1 of window 1
            return "Run pressed, long calculation confirmed"
        end try
        return "Run pressed"
    end tell
end tell
APPLESCRIPT
with_timeout 60 osascript "$OUT/run.applescript" || echo "run failed"
sleep 5
shot "e1b-running"
sleep 120
view_tab "e2-h3-summary" "Summary"
print_shot "e2-h3-summary"
view_tab "e3-h3-optimisation" "Optimisation"
print_shot "e3-h3-optimisation"
view_tab "e4-h3-frequencies" "Frequencies"
print_shot "e4-h3-frequencies"

# --- ΨLab 20.0: nuclear motion (the second kind of calculation) ---
echo "== 20.0: nuclear motion: the switch, ns1 = 20, np = 24, jmax = 0 (Lanczos), then jmax = 2"
view_tab "f0-electronic-before-switch" "Summary"
with_timeout 20 osascript "$OUT/view-menu.applescript" "Nuclear-Motion Calculation" || echo "could not switch to nuclear motion"
sleep 3
dump_tab "f1-nuclear-sidebar" ""
print_shot "f1-nuclear-sidebar"
with_timeout 60 osascript "$OUT/run.applescript" || echo "run failed"
sleep 30
view_tab "f2-nuclear-levels" "Energy Levels"
print_shot "f2-nuclear-levels"
view_tab "f3-nuclear-expectation" "Averaged Geometry"
print_shot "f3-nuclear-expectation"
view_tab "f3b-nuclear-matrix" "Matrix Structure"
print_shot "f3b-nuclear-matrix"
view_tab "f3c-nuclear-eigenvectors" "Eigenvectors"
print_shot "f3c-nuclear-eigenvectors"
view_tab "f4-nuclear-convergence" "Convergence"
print_shot "f4-nuclear-convergence"
# 21.0: the charts of the linear solves, below the first screen of the tab
scroll_shot "f4b-nuclear-convergence-linear" "Linear solves"
view_tab "f5-nuclear-output" "Nuclear-Motion Output"
# (the output file is an NSTextView, not in the text dump: the levels are checked instead)
grep -E "^4638\.3055$" "$OUT/f2-nuclear-levels.txt" > /dev/null && echo "nuclear levels: ZPE 4638.3055 cm-1 (ns1 = 20, np = 24) found" || echo "nuclear levels: ZPE 4638.3055 cm-1 NOT found"
with_timeout 20 osascript "$OUT/view-menu.applescript" "Electronic Structure Calculation" || echo "could not switch to electronic structure"
sleep 3
view_tab "f6-electronic-after-switch" "Summary"
print_shot "f6-electronic-after-switch"
pkill -x "ΨLab" || true
sleep 3
open -n "$APP" --args -appearanceChoice light -psilabCalculation nuclear -psilabNuclearNs1 10 -psilabNuclearNp 12 -psilabNuclearJmax 2 -psilabNuclearLevels 6 -psilabNuclearSolver jacobi-davidson -psilabNuclearAutoRun YES
sleep 25
view_tab "f7-nuclear-j2-levels" "Energy Levels"
print_shot "f7-nuclear-j2-levels"
view_tab "f8-nuclear-j2-expectation" "Averaged Geometry"
view_tab "f9-nuclear-j2-matrix" "Matrix Structure"
print_shot "f9-nuclear-j2-matrix"

# --- ΨLab 21.0: Radau coordinates, Rayleigh quotient iteration (MINRES), the Gaussian start vector
echo "== 21.0: nuclear motion in Radau coordinates, Rayleigh quotient iteration, Gaussian start vector"
pkill -x "ΨLab" || true
sleep 3
open -n "$APP" --args -appearanceChoice light -psilabCalculation nuclear -psilabNuclearNs1 16 -psilabNuclearNp 20 -psilabNuclearJmax 0 -psilabNuclearLevels 4 -psilabNuclearSolver rayleigh-quotient -psilabNuclearLinear minres -psilabNuclearCoordinates radau -psilabNuclearGaussian YES -psilabNuclearConvergenceScroll linear -psilabNuclearAutoRun YES
sleep 25
dump_tab "g1-nuclear-radau-sidebar" ""
print_shot "g1-nuclear-radau-sidebar"
view_tab "g2-nuclear-radau-levels" "Energy Levels"
print_shot "g2-nuclear-radau-levels"
grep -E "^4638\.3055$" "$OUT/g2-nuclear-radau-levels.txt" > /dev/null && echo "Radau levels: ZPE 4638.3055 cm-1 (ns1 = 16, np = 20) found" || echo "Radau levels: ZPE 4638.3055 cm-1 NOT found"
view_tab "g3-nuclear-radau-eigenvectors" "Eigenvectors"
print_shot "g3-nuclear-radau-eigenvectors"
# (-psilabNuclearConvergenceScroll linear: the tab opens at its charts of the linear solves)
view_tab "g4-nuclear-radau-convergence" "Convergence"
print_shot "g4-nuclear-radau-convergence"
view_tab "g5-nuclear-radau-geometry" "Averaged Geometry"
print_shot "g5-nuclear-radau-geometry"

echo "== relaunching $APP (HeH+ STO-3G, Szabo-Ostlund) for the sections of the earlier versions"
pkill -x "ΨLab" || true
sleep 3
open -n "$APP" --args -psilabMolecule HeH+ -psilabBasis STO-3G -psilabSzaboOstlund YES
sleep 15

layout "launch-heh"
echo "::group::accessibility tree"
with_timeout 40 osascript "$OUT/tab.applescript" "__TREE__" | tee "$OUT/tree.txt"
echo "::endgroup::"
dump_tab "01-summary" ""
view_tab "02-hartree-fock" "Hartree–Fock"
view_tab "03-full-ci" "Full-CI"
view_tab "04-potential-curve" "Potential Curve"

echo "== running the potential-curve scan"
echo "::group::accessibility tree (potential curve tab)"
with_timeout 40 osascript "$OUT/tab.applescript" "__TREE__" | tee "$OUT/tree-curve.txt"
echo "::endgroup::"
echo "::group::click Run scan"
with_timeout 90 osascript "$OUT/tab.applescript" "__SCAN__" 
echo "::endgroup::"
sleep 8
dump_tab "05-curve-after-scan" ""

echo "== menu commands"
echo "::group::menu items"
cat > "$OUT/menus.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        set out to ""
        repeat with m in {"File", "Edit"}
            try
                set names to name of every menu item of menu 1 of menu bar item (m as text) of menu bar 1
                set out to out & m & ": " & (names as text) & linefeed
            end try
        end repeat
        return out
    end tell
end tell
APPLESCRIPT
with_timeout 30 osascript "$OUT/menus.applescript" | tee "$OUT/menus.txt"
echo "::endgroup::"

echo "== Edit > Copy Summary as Text"
pbcopy < /dev/null
cat > "$OUT/copy.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        click menu item "Copy Summary as Text" of menu 1 of menu bar item "Edit" of menu bar 1
    end tell
end tell
APPLESCRIPT
with_timeout 30 osascript "$OUT/copy.applescript" || echo "could not click the menu item"
sleep 1
echo "::group::clipboard after Copy Summary as Text"
pbpaste | tee "$OUT/clipboard-summary.txt"
echo "::endgroup::"

echo "== File > Export Results as Text (accept the save panel with Return)"
cat > "$OUT/export-text.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        set frontmost to true
        click menu item "Export Results as Text…" of menu 1 of menu bar item "File" of menu bar 1
        delay 3
        key code 36
        delay 3
    end tell
end tell
APPLESCRIPT
with_timeout 40 osascript "$OUT/export-text.applescript" || echo "could not drive the export"
EXPORTED=$(find "$HOME" -name 'PsiLab-HeH+-STO3G-R1.4632.txt' -mmin -10 2>/dev/null | head -n 1)
echo "::group::exported text file"
if [ -n "$EXPORTED" ]; then
    echo "file: $EXPORTED ($(wc -l < "$EXPORTED") lines)"
    cp "$EXPORTED" "$OUT/exported-results.txt"
    head -n 40 "$EXPORTED"
    echo "[...]"
    tail -n 60 "$EXPORTED"
else
    echo "exported file NOT found"
fi
echo "::endgroup::"

# --- dark mode pass: relaunch with the in-app appearance forced to dark ---------------
echo "== dark mode pass"
pkill -x "ΨLab" || true
sleep 3
open -n "$APP" --args -appearanceChoice dark -psilabMolecule HeH+ -psilabBasis STO-3G -psilabSzaboOstlund YES
sleep 15
click_shot() {   # $1 = file prefix, $2 = tab to click ("" = none)
    with_timeout 60 osascript "$OUT/tab.applescript" "$2" > "$OUT/$1.txt" 2>&1
    shot "$1"
}
view_shot() {   # $1 = file prefix, $2 = View menu item
    with_timeout 20 osascript "$OUT/view-menu.applescript" "$2" > /dev/null 2>&1
    sleep 2
    shot "$1"
}
view_shot "dark-01-summary" "Summary"
view_shot "dark-02-hartree-fock" "Hartree–Fock"
view_shot "dark-03-full-ci" "Full-CI"
view_shot "dark-04-potential-curve" "Potential Curve"
with_timeout 90 osascript "$OUT/tab.applescript" "__SCAN__" | grep -E "CLICKED|NOT FOUND"
sleep 8
shot "dark-05-curve-after-scan"

# --- ΨLab 3.0: H2, LiH, H2O (launch arguments select molecule, basis set and method) ---------
relaunch() {   # $@ = launch arguments
    pkill -x "ΨLab" || true
    sleep 3
    open -n "$APP" --args "$@"
    sleep 20
}
echo "== H2O, RHF/aug-cc-pVDZ"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF
view_tab "10-h2o-rhf-summary" "Summary"
layout "h2o-rhf-relaunch"
view_tab "11-h2o-rhf-hartree-fock" "Hartree–Fock"
view_tab "12-h2o-rhf-output" "Output"

echo "== LiH, Full CI/aug-cc-pVDZ"
relaunch -appearanceChoice light -psilabMolecule LiH -psilabBasis aug-cc-pVDZ -psilabMethod FCI
view_tab "13-lih-fci-summary" "Summary"
view_tab "14-lih-fci-full-ci" "Full-CI"

echo "== H2, Full CI/aug-cc-pVDZ potential curve"
relaunch -appearanceChoice light -psilabMolecule H2 -psilabBasis aug-cc-pVDZ -psilabMethod FCI
view_tab "15-h2-curve" "Potential Curve"
with_timeout 90 osascript "$OUT/tab.applescript" "__SCAN__" | grep -E "CLICKED|NOT FOUND"
sleep 25
dump_tab "16-h2-curve-after-scan" ""

echo "== H2O, Full CI (must be refused) and a long calculation (warning)"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod FCI
view_tab "17-h2o-fci-refused" "Summary"
relaunch -appearanceChoice dark -psilabMolecule H2O -psilabBasis aug-cc-pCV5Z -psilabMethod RHF
view_tab "18-h2o-cv5z-warning-dark" "Summary"

# --- ΨLab 4.0: coupled cluster ------------------------------------------------------------
cat > "$OUT/method-menu.applescript" <<'APPLESCRIPT'
on run argv
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            set els to entire contents of window 1
            repeat with e in els
                try
                    if (role of e) is "AXPopUpButton" and (name of e) is "Method" then
                        click e
                        return "CLICKED Method pop-up"
                    end if
                end try
            end repeat
        end tell
    end tell
    return "NOT FOUND Method pop-up"
end run
APPLESCRIPT
echo "== H2O, CCSD(T)/aug-cc-pVDZ and the method menu"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod "CCSD(T)"
view_tab "19-h2o-ccsd-t-summary" "Summary"
view_tab "20-h2o-ccsd-t-coupled-cluster" "Coupled Cluster"
with_timeout 30 osascript "$OUT/method-menu.applescript" | grep -E "CLICKED|NOT FOUND"
sleep 2
shot "21-method-menu"
osascript -e 'tell application "System Events" to key code 53' || true

echo "== LiH, CCSDT/aug-cc-pVDZ (dark) and H2O CCSDTQ/aug-cc-pVTZ (refused)"
relaunch -appearanceChoice dark -psilabMolecule LiH -psilabBasis aug-cc-pVDZ -psilabMethod CCSDT
view_tab "22-lih-ccsdt-coupled-cluster-dark" "Coupled Cluster"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVTZ -psilabMethod CCSDTQ
view_tab "23-h2o-ccsdtq-refused" "Summary"


# --- ΨLab 5.0: job types ------------------------------------------------------------------
cat > "$OUT/job-menu.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        set frontmost to true
        set els to entire contents of window 1
        repeat with e in els
            try
                set r to role of e
                set n to ""
                set d to ""
                try
                    set n to name of e
                end try
                try
                    set d to description of e
                end try
                if (r is "AXMenuButton" or r is "AXPopUpButton") and (n is "Job" or d is "Job") then
                    click e
                    return "CLICKED Job menu (" & r & ")"
                end if
            end try
        end repeat
    end tell
end tell
return "NOT FOUND Job menu"
APPLESCRIPT
cat > "$OUT/cancel.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        set frontmost to true
        keystroke "." using command down
    end tell
end tell
return "SENT Cmd-."
APPLESCRIPT
# Only the Output tab's text view holds the output file: print its first lines from the tree.
cat > "$OUT/output-text.applescript" <<'APPLESCRIPT'
tell application "System Events"
    tell process "ΨLab"
        set els to entire contents of window 1
        repeat with e in els
            try
                if (role of e) is "AXTextArea" then
                    set v to value of e as text
                    return "OUTPUT FILE: " & (length of v) & " characters" & linefeed & (text 1 thru (my min_len(v)) of v)
                end if
            end try
        end repeat
    end tell
end tell
return "NO TEXT AREA"
on min_len(v)
    if (length of v) < 6000 then return length of v
    return 6000
end min_len
APPLESCRIPT

echo "== job menu (H2O, RHF)"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob opt
with_timeout 30 osascript "$OUT/job-menu.applescript"
sleep 2
shot "30-job-menu"
osascript -e 'tell application "System Events" to key code 53' || true

echo "== H2O RHF/aug-cc-pVDZ geometry optimisation from R = 1.0 Å, 110°"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob opt -psilabR 1.0 -psilabAngle 110
sleep 10
view_tab "31-h2o-rhf-opt-optimisation" "Optimisation"
layout "h2o-optimisation"
# ΨLab 6.0: the molecule of the shown step (the last); still until Play steps is pressed
frames "31-opt-still" 1.8
print_shot "31-h2o-rhf-opt-optimisation"
view_tab "32-h2o-rhf-opt-summary" "Summary"
# Play steps, switched on by a launch argument (a scripted click on the button takes longer than
# the AppleScript time limit on this page): the shown step changes while it plays.
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob opt -psilabR 1.0 -psilabAngle 110 -psilabPlaySteps YES
sleep 10
view_tab "31b-h2o-opt-playing" "Optimisation"
frames "31-opt-playing" 1.8
dump_tab "31c-h2o-opt-playing-later" ""
print_shot "31b-h2o-opt-playing"

echo "== LiH CCSD(T)/aug-cc-pVDZ frequencies and thermochemistry"
relaunch -appearanceChoice light -psilabMolecule LiH -psilabBasis aug-cc-pVDZ -psilabMethod "CCSD(T)" -psilabJob freq
sleep 15
view_tab "33-lih-ccsd-t-frequencies" "Frequencies"
frames "33-lih-mode"
print_shot "33-lih-ccsd-t-frequencies"
view_tab "34-lih-ccsd-t-thermochemistry" "Thermochemistry"

echo "== H2O opt+freq+SP: RHF/aug-cc-pVDZ geometry and frequencies, CCSD(T)/aug-cc-pVTZ single point"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob opt+freq+sp -psilabSPMethod "CCSD(T)" -psilabSPBasis aug-cc-pVTZ -psilabR 1.0 -psilabAngle 110 -psilabAutoRun YES
sleep 60
view_tab "35-h2o-job5-summary" "Summary"
view_tab "35b-h2o-job5-frequencies" "Frequencies"
frames "35b-h2o-modes"
print_shot "35b-h2o-job5-frequencies"
view_tab "36-h2o-job5-thermochemistry" "Thermochemistry"
view_tab "37-h2o-job5-output" "Output"
echo "::group::output file (Output tab)"
with_timeout 60 osascript "$OUT/output-text.applescript" | head -c 7000
echo
echo "::endgroup::"

echo "== Output tab while a job runs (H2O CCSD(T)/aug-cc-pVTZ opt+freq), then Cancel"
relaunch -appearanceChoice dark -psilabMolecule H2O -psilabBasis aug-cc-pVTZ -psilabMethod "CCSD(T)" -psilabJob opt+freq -psilabAutoRun YES -psilabTab output
sleep 10
shot "38-output-running-dark"
sleep 15
shot "39-output-running-later-dark"
with_timeout 20 osascript "$OUT/cancel.applescript"
sleep 5
dump_tab "40-output-cancelled-dark" ""

echo "== H2O CCSDTQ frequencies (refused)"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod CCSDTQ -psilabJob freq
view_tab "41-h2o-ccsdtq-freq-refused" "Summary"

# --- ΨLab 7.0: HeH+ in a correlation-consistent basis set, H2O2 and NH3, relaxed scans and
# transition-state searches -------------------------------------------------------------------
echo "== HeH+ Full CI/aug-cc-pVTZ opt+freq (the general path)"
relaunch -appearanceChoice light -psilabMolecule HeH+ -psilabBasis aug-cc-pVTZ -psilabMethod FCI -psilabJob opt+freq -psilabAutoRun YES
sleep 15
view_tab "50-heh-fci-opt-freq-summary" "Summary"
view_tab "51-heh-fci-frequencies" "Frequencies"

echo "== H2O2 RHF/aug-cc-pVDZ relaxed scan of the H-O-O-H torsion, 0 ... 180 degrees in 30-degree steps"
# one launch: the points are played (Play points on), so the frames differ
relaunch -appearanceChoice light -psilabMolecule H2O2 -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob scan -psilabScan "D 3 1 2 4 S 6 30" -psilabScanStart 0 -psilabAutoRun YES -psilabTab scan -psilabPlaySteps YES
shot "52-h2o2-scan-running"
sleep 80
view_tab "53-h2o2-relaxed-scan" "Relaxed Scan"
frames "53-scan-playing" 1.8
print_shot "53-h2o2-relaxed-scan"
view_tab "54-h2o2-scan-summary" "Summary"
view_tab "55-h2o2-scan-output" "Output"

echo "== NH3 RHF/aug-cc-pVDZ inversion transition state (ts+freq from H-N-H 119 degrees)"
relaunch -appearanceChoice light -psilabMolecule NH3 -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob ts+freq -psilabAngle 119 -psilabAutoRun YES
sleep 40
view_tab "57-nh3-ts-optimisation" "Optimisation"
print_shot "57-nh3-ts-optimisation"
view_tab "58-nh3-ts-frequencies" "Frequencies"
frames "58-nh3-imaginary-mode"
print_shot "58-nh3-ts-frequencies"
view_tab "59-nh3-ts-summary" "Summary"

echo "== H2O2 RHF/aug-cc-pVDZ trans transition state with CalcAll (from 160 degrees), sidebar of the TS job"
relaunch -appearanceChoice light -psilabMolecule H2O2 -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob ts+freq -psilabDihedral 160 -psilabTSHessian calcall -psilabAutoRun YES
sleep 60
view_tab "60-h2o2-ts-calcall-optimisation" "Optimisation"
view_tab "61-h2o2-ts-frequencies" "Frequencies"
print_shot "61-h2o2-ts-frequencies"

# --- ΨLab 8.0: the molecular orbitals ------------------------------------------------------------
echo "== orbitals: H2O RHF/aug-cc-pVDZ HOMO and LUMO isosurfaces, NH3 contour plot, HeH+ Full CI natural orbitals, Szabo-Ostlund HeH+"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabTab orbitals
sleep 8
dump_tab "64-h2o-orbitals-homo" ""
print_shot "64-h2o-orbitals-homo"
relaunch -appearanceChoice dark -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabTab orbitals -psilabOrbital LUMO -psilabIsovalue 0.02
sleep 8
shot "65-h2o-orbitals-lumo-dark"
relaunch -appearanceChoice light -psilabMolecule NH3 -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabTab orbitals -psilabOrbitalView contour
sleep 8
dump_tab "66-nh3-orbital-contour" ""
print_shot "66-nh3-orbital-contour"
relaunch -appearanceChoice light -psilabMolecule HeH+ -psilabBasis aug-cc-pVDZ -psilabMethod FCI -psilabTab orbitals -psilabOrbitalSet natural
sleep 8
dump_tab "67-heh-fci-natural-orbitals" ""
relaunch -appearanceChoice light -psilabMolecule HeH+ -psilabBasis STO-3G -psilabSzaboOstlund YES -psilabTab orbitals
sleep 5
dump_tab "68-heh-sto3g-orbitals" ""
view_tab "69-orbitals-view-menu" "Orbitals"

# --- ΨLab 9.0: custom molecules ------------------------------------------------------------------
# The sandboxed app reads files only through the file panel, so the structures are passed as text
# (-psilabStructure); the molecule is kept in the molecule menu by its name (CO$_2$ drawn as CO₂).
cat > "$OUT/file-menu.applescript" <<'APPLESCRIPT'
on run argv
    tell application "System Events"
        tell process "ΨLab"
            set frontmost to true
            click menu item (item 1 of argv) of menu 1 of menu bar item "File" of menu bar 1
        end tell
    end tell
    return "FILE MENU " & (item 1 of argv)
end run
APPLESCRIPT
echo "== custom molecule: CO2 from an XYZ text (name CO\$_2\$), RHF/aug-cc-pVDZ opt+freq, optimised geometry kept"
CO2_XYZ=$'3\nCO$_2$\nC   0.000000   0.000000   0.000000\nO   0.000000   0.000000   1.170000\nO   0.000000   0.000000  -1.170000\n'
relaunch -appearanceChoice light -psilabStructure "$CO2_XYZ" -psilabStructureName co2.xyz -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob opt+freq -psilabAutoRun YES -psilabSaveOptimized YES
sleep 40
view_tab "70-co2-custom-summary" "Summary"
print_shot "70-co2-custom-summary"
view_tab "71-co2-custom-optimisation" "Optimisation"
print_shot "71-co2-custom-optimisation"
view_tab "72-co2-custom-frequencies" "Frequencies"
frames "72-co2-bending-mode"
view_tab "73-co2-custom-output" "Output"

echo "== custom molecule: methane from a Gaussian input (charge 0, multiplicity 1), RHF frequencies; the menu lists CO2 and its optimised geometry"
CH4_GJF=$'%chk=ch4.chk\n#p HF/aug-cc-pVDZ freq\n\nmethane, CH$_4$\n\n0 1\nC 0 0 0\nH 0.629118 0.629118 0.629118\nH -0.629118 -0.629118 0.629118\nH -0.629118 0.629118 -0.629118\nH 0.629118 -0.629118 -0.629118\n\n'
relaunch -appearanceChoice dark -psilabStructure "$CH4_GJF" -psilabStructureName ch4.gjf -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob freq -psilabAutoRun YES
sleep 30
view_tab "74-ch4-custom-summary" "Summary"
view_tab "75-ch4-custom-frequencies" "Frequencies"
print_shot "75-ch4-custom-frequencies"
with_timeout 20 osascript "$OUT/tab.applescript" "Atom/Molecule" > /dev/null 2>&1 || true
sleep 2
shot "76-molecule-menu-with-custom-molecules"
print_shot "76-molecule-menu-with-custom-molecules"
osascript -e 'tell application "System Events" to key code 53' || true
echo "== File > Open Structure File… (the file panel; closed with Escape)"
with_timeout 20 osascript "$OUT/file-menu.applescript" "Open Structure File…" || echo "could not open the file panel"
sleep 3
shot "77-open-structure-panel"
osascript -e 'tell application "System Events" to key code 53' || true
sleep 1

echo "== custom molecule refused: an odd number of electrons (OH radical) and a heavy element (RbH)"
relaunch -appearanceChoice light -psilabStructure $'2\nOH\nO 0 0 0\nH 0 0 0.97\n' -psilabStructureName oh.xyz -psilabBasis aug-cc-pVDZ -psilabMethod RHF
dump_tab "78-oh-radical-refused" ""
relaunch -appearanceChoice light -psilabStructure $'2\nRbH\nRb 0 0 0\nH 0 0 2.37\n' -psilabStructureName rbh.xyz -psilabBasis aug-cc-pVDZ -psilabMethod RHF
dump_tab "79-rbh-refused" ""

# --- ΨLab 10.0: analytic Hessians of Full CI, CCSD and CCSD(T) -----------------------------------
echo "== H2O CCSD/aug-cc-pVDZ frequencies with the analytic Hessian (sidebar: the Hessian section)"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod CCSD -psilabJob freq -psilabAutoRun YES
dump_tab "80-h2o-ccsd-analytic-hessian-sidebar" ""
sleep 60
view_tab "81-h2o-ccsd-analytic-frequencies" "Frequencies"
print_shot "81-h2o-ccsd-analytic-frequencies"
view_tab "82-h2o-ccsd-analytic-output" "Output"
echo "== the same with the numerical Hessian (estimate only)"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod "CCSD(T)" -psilabJob freq -psilabHessian numerical
dump_tab "83-h2o-ccsd-t-numerical-hessian-sidebar" ""

echo "== NH3 scan settings in the sidebar (umbrella) and the refused TS search of a diatomic"
relaunch -appearanceChoice light -psilabMolecule NH3 -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob scan
dump_tab "62-nh3-scan-sidebar" ""
relaunch -appearanceChoice light -psilabMolecule LiH -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob ts
view_tab "63-lih-ts-refused" "Summary"

# --- ΨLab 12.0: the 3D molecule viewer -------------------------------------------------------------
# Its display settings are stored defaults (viewer3D.*), so launch arguments set them: sticks with
# numbered atom labels for the optimisation of H2O2, van der Waals spheres on black for the
# animated normal modes of NH3 (two captures must differ), the Summary of H2O2 in the default style.
echo "== 3D viewer: H2O2 RHF optimisation in sticks with atom numbers, NH3 modes as spheres on black"
relaunch -appearanceChoice light -psilabMolecule H2O2 -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob opt -psilabAutoRun YES -viewer3D.representation sticks -viewer3D.labels numbered
sleep 25
view_tab "80-h2o2-opt-3d-sticks" "Optimisation"
print_shot "80-h2o2-opt-3d-sticks"
view_tab "81-h2o2-summary-3d" "Summary"
print_shot "81-h2o2-summary-3d"
relaunch -appearanceChoice dark -psilabMolecule NH3 -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob freq -psilabAutoRun YES -viewer3D.representation spheres -viewer3D.background black
sleep 25
view_tab "82-nh3-freq-3d-spheres-dark" "Frequencies"
frames "82-nh3-spheres-mode"
print_shot "82-nh3-freq-3d-spheres-dark"

# --- ΨLab 13.0: an output file opened again ----------------------------------------------------------
# examples/outputs/h2o-ccsd-opt-freq.out, written by psilab 13.0: its archive section restores every
# tab without a calculation (the orbitals are computed again from the basis set and coefficients).
echo "== open an output file of psilab (H2O CCSD/aug-cc-pVDZ opt+freq): Summary, Optimisation, Frequencies, Orbitals"
# The app is sandboxed: a path given on the command line must lie in its container (a file chosen
# in File > Open Output File… is granted by the file panel). The container exists after the
# launches above.
CONTAINER="$HOME/Library/Containers/org.psilab.PsiLab/Data"
mkdir -p "$CONTAINER/tmp"
cp examples/outputs/h2o-ccsd-opt-freq.out "$CONTAINER/tmp/" || echo "could not copy the example output into the container"
EXAMPLE_OUT="$CONTAINER/tmp/h2o-ccsd-opt-freq.out"
relaunch -appearanceChoice light -psilabOpenOutput "$EXAMPLE_OUT"
view_tab "90-opened-output-summary" "Summary"
print_shot "90-opened-output-summary"
view_tab "91-opened-output-optimisation" "Optimisation"
view_tab "92-opened-output-frequencies" "Frequencies"
print_shot "92-opened-output-frequencies"
relaunch -appearanceChoice light -psilabOpenOutput "$EXAMPLE_OUT" -psilabTab orbitals
sleep 5
dump_tab "93-opened-output-orbitals" ""
print_shot "93-opened-output-orbitals"

# --- ΨLab 22.0: a nuclear-motion output file opened again (its archive: every nuclear-motion tab) -----
# examples/outputs/water-radau-small.out, written by psilab 22.0 from examples/water-radau-small.inp
echo "== 22.0: open a nuclear-motion output file of psilab (H2(16)O, Radau, ns1 = 10, np = 12, 12 levels)"
cp examples/outputs/water-radau-small.out "$CONTAINER/tmp/" || echo "could not copy the nuclear-motion output into the container"
relaunch -appearanceChoice light -psilabOpenOutput "$CONTAINER/tmp/water-radau-small.out"
sleep 3
view_tab "94-opened-nuclear-levels" "Energy Levels"
print_shot "94-opened-nuclear-levels"
grep -E "^4640\.5477$" "$OUT/94-opened-nuclear-levels.txt" > /dev/null && echo "opened nuclear levels: ZPE 4640.5477 cm-1 found" || echo "opened nuclear levels: ZPE 4640.5477 cm-1 NOT found"
view_tab "95-opened-nuclear-geometry" "Averaged Geometry"
print_shot "95-opened-nuclear-geometry"
view_tab "96-opened-nuclear-matrix" "Matrix Structure"
print_shot "96-opened-nuclear-matrix"
view_tab "97-opened-nuclear-eigenvectors" "Eigenvectors"
print_shot "97-opened-nuclear-eigenvectors"
view_tab "98-opened-nuclear-convergence" "Convergence"
view_tab "99-opened-nuclear-output" "Nuclear-Motion Output"

# --- ΨLab 14.0: the Pople 6-31G and 6-311G families and the def2 basis sets --------------------------
# 6-31G* has Cartesian d functions (6D): its orbitals pass through the app with "pure": false.
echo "== 14.0: H2O RHF/6-31G* orbitals (Cartesian d), NH3 RHF/def2-TZVP opt+freq, H2O CCSD(T)/6-311++G(3df,3pd)"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis "6-31G*" -psilabMethod RHF -psilabTab orbitals
sleep 8
dump_tab "95-h2o-631gs-orbitals" ""
print_shot "95-h2o-631gs-orbitals"
relaunch -appearanceChoice light -psilabMolecule NH3 -psilabBasis def2-TZVP -psilabMethod RHF -psilabJob opt+freq -psilabAutoRun YES
sleep 25
view_tab "96-nh3-def2-tzvp-summary" "Summary"
print_shot "96-nh3-def2-tzvp-summary"
view_tab "97-nh3-def2-tzvp-frequencies" "Frequencies"
relaunch -appearanceChoice dark -psilabMolecule H2O -psilabBasis "6-311++G(3df,3pd)" -psilabMethod "CCSD(T)" -psilabAutoRun YES -psilabTab output
sleep 20
view_tab "98-h2o-6-311ppg-3df-3pd-output" "Output"

# Elements beyond Kr: the def2 sets with their effective core potentials (HI, ECP28 on I; PbO, ECP60)
echo "== 14.0: HI RHF/def2-TZVP opt+freq (I with the def2 ECP), PbO CCSD(T)/def2-SVP energy"
relaunch -appearanceChoice light -psilabStructure $'2\nHI\nH 0 0 0\nI 0 0 1.62\n' -psilabStructureName hi.xyz -psilabBasis def2-TZVP -psilabMethod RHF -psilabJob opt+freq -psilabAutoRun YES
sleep 20
view_tab "99-hi-def2-tzvp-summary" "Summary"
print_shot "99-hi-def2-tzvp-summary"
view_tab "99b-hi-def2-tzvp-frequencies" "Frequencies"
relaunch -appearanceChoice dark -psilabStructure $'2\nPbO\nPb 0 0 0\nO 0 0 1.92\n' -psilabStructureName pbo.xyz -psilabBasis def2-SVP -psilabMethod "CCSD(T)" -psilabAutoRun YES -psilabTab output
sleep 20
view_tab "99c-pbo-ccsdt-output" "Output"

# --- ΨLab 15.0: open shells, the atoms and new molecules, MPn/QCI, optimiser options -------------
echo "== 15.0: O2 UMP2/6-31G* opt+freq (RecalcFC=2, GEDIIS), N atom ROHF, H2O QCISD(T)/6-31G opt+freq, HeH+ STO-3G"
relaunch -appearanceChoice light -psilabMolecule O2 -psilabBasis "6-31G*" -psilabMethod MP2 -psilabJob opt+freq -psilabOptHessian recalcfc -psilabRecalcFC 2 -psilabGEDIIS YES -psilabAutoRun YES
sleep 25
view_tab "a0-o2-ump2-summary" "Summary"
print_shot "a0-o2-ump2-summary"
view_tab "a1-o2-ump2-optimisation" "Optimisation"
print_shot "a1-o2-ump2-optimisation"
view_tab "a1b-o2-ump2-frequencies" "Frequencies"
print_shot "a1b-o2-ump2-frequencies"
# 16.0: the segmented controls of the sidebar (Hessian: Initial | CalcFC | CalcAll | RecalcFC) lie
# inside the 340-pt panel
layout "a1d-sidebar-controls"
grep -E "Initial|CalcFC|RecalcFC" "$OUT/layout-a1d-sidebar-controls.txt" || echo "Hessian control not in the accessibility layout"
# the molecule menu of 15.0: the molecules (built in and with fixed geometries) and the atoms
with_timeout 20 osascript "$OUT/tab.applescript" "Atom/Molecule" > /dev/null 2>&1 || true
sleep 2
shot "a1c-molecule-menu-15"
print_shot "a1c-molecule-menu-15"
osascript -e 'tell application "System Events" to key code 53' || true
sleep 1
relaunch -appearanceChoice dark -psilabMolecule N -psilabBasis "6-31G*" -psilabMethod HF -psilabReference ROHF -psilabAutoRun YES -psilabTab output
sleep 10
view_tab "a2-n-rohf-output" "Output"
relaunch -appearanceChoice light -psilabMolecule H2O -psilabBasis 6-31G -psilabMethod "QCISD(T)" -psilabJob opt+freq -psilabAutoRun YES
sleep 40
view_tab "a3-h2o-qcisdt-frequencies" "Frequencies"
print_shot "a3-h2o-qcisdt-frequencies"
relaunch -appearanceChoice light -psilabMolecule HeH+ -psilabBasis STO-3G -psilabSzaboOstlund YES
sleep 5
view_tab "a4-heh-sto3g-summary" "Summary"
print_shot "a4-heh-sto3g-summary"
# H2O2 relaxed scan of the torsion (the larger 3D viewer of a scan point) and NO2 ROMP2 opt+freq
relaunch -appearanceChoice light -psilabMolecule H2O2 -psilabBasis 6-31G -psilabMethod HF -psilabJob scan -psilabScan "D 3 1 2 4 S 6 30" -psilabScanStart 0 -psilabAutoRun YES
sleep 30
view_tab "a5-h2o2-scan-15" "Relaxed Scan"
print_shot "a5-h2o2-scan-15"
relaunch -appearanceChoice dark -psilabMolecule NO2 -psilabBasis STO-3G -psilabMethod MP2 -psilabReference ROHF -psilabJob opt+freq -psilabAutoRun YES
sleep 40
view_tab "a6-no2-romp2-summary" "Summary"
print_shot "a6-no2-romp2-summary"
view_tab "a7-no2-romp2-frequencies" "Frequencies"
print_shot "a7-no2-romp2-frequencies"
# 16.0: the estimate of an impossible job (Full CI of Cl2 in aug-cc-pCV5Z) shows its time and memory
relaunch -appearanceChoice light -psilabMolecule Cl2 -psilabBasis aug-cc-pCV5Z -psilabMethod FCI
dump_tab "a8-cl2-fci-impossible-estimate" ""
grep -E "Estimated time|Full-CI memory|Determinants" "$OUT/a8-cl2-fci-impossible-estimate.txt" || echo "estimate rows not found"

# --- ΨLab 17.0: the start-up example (no launch arguments), the Summary viewer with its caption,
# the save buttons of the optimisation steps, and HeH+ STO-3G with the Szabo-Ostlund switch off ---
echo "== 17.0: launch without arguments: H2 Full CI/aug-cc-pCV5Z//Full CI/aug-cc-pCVQZ (opt+freq+sp)"
pkill -x "ΨLab" || true
sleep 3
open -n "$APP" --args -appearanceChoice light
sleep 15
view_tab "b1-startup-example-summary" "Summary"
print_shot "b1-startup-example-summary"
grep -E "//|Drag to turn" "$OUT/b1-startup-example-summary.txt" || echo "Level1//Level2 or caption not found"
# 18.0: the start-up example has no banner (the buttons are in the accessibility tree)
with_timeout 60 osascript "$OUT/tab.applescript" "__TREE__" > "$OUT/b1-startup-example-tree.txt" 2>&1
if grep -E "Example: H|Results read from" "$OUT/b1-startup-example-tree.txt"; then echo "unexpected example banner"; else echo "no example banner (expected)"; fi
view_tab "b2-startup-example-optimisation" "Optimisation"
print_shot "b2-startup-example-optimisation"
with_timeout 60 osascript "$OUT/tab.applescript" "__TREE__" > "$OUT/b2-startup-example-tree.txt" 2>&1
grep -E "Save (Step|All Steps) as|this step as|all steps as" "$OUT/b2-startup-example-tree.txt" || echo "save buttons not found"
view_tab "b3-startup-example-frequencies" "Frequencies"
relaunch -appearanceChoice light -psilabMolecule HeH+ -psilabBasis STO-3G
sleep 5
dump_tab "b4-heh-sto3g-szabo-ostlund-off" ""
grep -E "Szabo–Ostlund calculation" "$OUT/b4-heh-sto3g-szabo-ostlund-off.txt" || echo "Szabo-Ostlund switch not found"

# --- ΨLab 18.0: the complete starting geometry of the molecules of the menu ---
# (the derived rows are static text: they are in the text dump, not in the accessibility tree)
echo "== 18.0: starting geometry of H3+, CH4, N2O and N2"
for spec in "H3+|R\(H1–H2\), R\(H1–H3\) \(Å\)|Angle H2–H1–H3 \(°\)|R\(H2–H3\) \(Å\)" \
            "CH4|R\(C–H\) \(Å\)|Angle H–C–H \(°\)" \
            "N2O|R\(N–N\) \(Å\)|R\(N–O\) \(Å\)|Angle N–N–O \(°\)" \
            "N2|R\(N–N\) \(Å\)"; do
    mol="${spec%%|*}"
    relaunch -appearanceChoice light -psilabMolecule "$mol" -psilabJob opt
    sleep 3
    dump_tab "c1-geometry-$mol" ""
    rest="${spec#*|}"
    IFS='|' read -r -a rows <<< "$rest"
    for row in "${rows[@]}"; do
        grep -E "^$row" "$OUT/c1-geometry-$mol.txt" > /dev/null && echo "$mol: found $row" || echo "$mol: geometry row $row not found"
    done
done

# --- ΨLab 19.0: H3+ in the 3D viewers (three bonds, H1 at the apex), atom labels as badges ---
echo "== 19.0: H3+ Full CI/aug-cc-pVDZ opt+freq+sp, numbered atom labels"
relaunch -appearanceChoice dark -viewer3D.labels numbered -psilabMolecule H3+ -psilabBasis aug-cc-pVDZ -psilabMethod FCI -psilabJob opt+freq+sp -psilabAutoRun YES
sleep 40
view_tab "d1-h3-optimisation" "Optimisation"
print_shot "d1-h3-optimisation"
view_tab "d2-h3-frequencies" "Frequencies"
print_shot "d2-h3-frequencies"
grep -E "mdyn/Å \(1 mdyn/Å = 100 N/m|974.88 km/mol" "$OUT/d2-h3-frequencies.txt" || echo "unit note not found"
relaunch -appearanceChoice light -viewer3D.labels numbered -psilabMolecule H2O -psilabBasis aug-cc-pVDZ -psilabMethod RHF -psilabJob opt -psilabR 1.0 -psilabAngle 110 -psilabAutoRun YES
sleep 10
view_tab "d3-h2o-optimisation-labels-light" "Optimisation"
print_shot "d3-h2o-optimisation-labels-light"

# --- ΨLab 23.0: the nuclear-motion startup example (H2(16)O, J = 10, ns1 = np = 50, 100 levels) in
# every tab, the one-line headers of the Geometry tab, the square R1-theta map, and a new tab with
# its own calculation (File > New Tab) ---
echo "== 23.0: launch without arguments, the nuclear-motion interface: the startup example"
relaunch -appearanceChoice light
sleep 10
view_tab "h1-nuclear-example-levels" "Nuclear-Motion Calculation"
print_shot "h1-nuclear-example-levels"
grep -E "5931\.37|J = 10" "$OUT/h1-nuclear-example-levels.txt" | head -3 || true
grep -E "5931\.37" "$OUT/h1-nuclear-example-levels.txt" > /dev/null && echo "nuclear example: E1 5931.37 cm-1 found" || echo "nuclear example: E1 5931.37 cm-1 NOT found"
view_tab "h2-nuclear-example-geometry" "Averaged Geometry"
print_shot "h2-nuclear-example-geometry"
view_tab "h3-nuclear-example-matrix" "Matrix Structure"
view_tab "h4-nuclear-example-eigenvectors" "Eigenvectors"
print_shot "h4-nuclear-example-eigenvectors"
dump_tab "h5-nuclear-example-r1-theta" "R₁ – θ"
print_shot "h5-nuclear-example-r1-theta"
view_tab "h6-nuclear-example-convergence" "Convergence"
view_tab "h7-nuclear-example-output" "Nuclear-Motion Output"
echo "== 23.0: File > New Tab: a second calculation (the electronic-structure example) in a tab"
with_timeout 20 osascript "$OUT/file-menu.applescript" "New Tab" || echo "could not open a new tab"
sleep 10
dump_tab "h8-new-tab" ""
print_shot "h8-new-tab"
grep -E "Electronic structure" "$OUT/h8-new-tab.txt" > /dev/null && echo "new tab: electronic-structure interface" || echo "new tab: electronic-structure text NOT found"

# 20.0/21.0: the screenshots of the 20.0 and 21.0 sections once more at the end of the log (the job
# log read through the API is its last part only; the sections run first)
echo "== 20.0/21.0 screenshots again"
for name in e1-h3-chosen e2-h3-summary e3-h3-optimisation e4-h3-frequencies f1-nuclear-sidebar f2-nuclear-levels \
            f3-nuclear-expectation f3b-nuclear-matrix f3c-nuclear-eigenvectors f4-nuclear-convergence \
            f4b-nuclear-convergence-linear \
            f6-electronic-after-switch f7-nuclear-j2-levels f9-nuclear-j2-matrix g1-nuclear-radau-sidebar \
            g2-nuclear-radau-levels g3-nuclear-radau-eigenvectors g4-nuclear-radau-convergence g5-nuclear-radau-geometry; do
    [ -f "$OUT/$name.png" ] && print_shot "$name" || echo "no screenshot $name"
done
grep -E "^4638\.3055$" "$OUT/f2-nuclear-levels.txt" > /dev/null && echo "nuclear levels: ZPE 4638.3055 cm-1 found" || echo "nuclear levels: ZPE 4638.3055 cm-1 NOT found"

pkill -x "ΨLab" || true

# Crash reports of this run (none expected): print the crashed thread's backtrace.
echo "::group::crash reports"
find "$HOME/Library/Logs/DiagnosticReports" /Library/Logs/DiagnosticReports -name '*Lab*' -mmin -60 2>/dev/null | while read -r report; do
    echo "== $report"
    python3 - "$report" <<'PY' || head -c 6000 "$report"
import json, sys
text = open(sys.argv[1]).read()
header, _, body = text.partition("\n")
try:
    d = json.loads(body)
except Exception:
    print(text[:6000]); sys.exit(0)
print("exception:", d.get("exception"), "termination:", d.get("termination", {}).get("indicator"))
print("asi:", d.get("asi"))
images = d.get("usedImages", [])
for t in d.get("threads", []):
    if t.get("triggered"):
        for f in t.get("frames", [])[:40]:
            img = images[f.get("imageIndex", 0)].get("name", "?") if images else "?"
            print("  ", img, f.get("symbol", "?"), "+", f.get("symbolLocation", ""))
PY
done
echo "::endgroup::"
exit 0
