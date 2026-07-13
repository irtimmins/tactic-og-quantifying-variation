*==============================================================================
* profile_cohort.do
*
* Describes any health dataset. For each column it works out
* what it holds, then prints a summary that suits it.
* It never prints individual rows. For columns that hold patient identifiers
* or long notes it reports only how many there are and how long they are,
* not what they say.
*==============================================================================

version 14.0
clear all
set more off

* Which file to read, and where to save the write-up. Passing a filename when
* you run the script overrides the first of these.
local infile `"`1'"'
if `"`infile'"' == "" ///
    local infile "D:/Projects/#2045_ICON_TACTIC/Project4_OG_variation_deviants/tactic-og-quantifying-variation/Data/test_data/og_ncras_test.dta"

if `"$logfile"' == "" ///
    global logfile "D:/Projects/#2045_ICON_TACTIC/Project4_OG_variation_deviants/tactic-og-quantifying-variation/Data/test_data/profile_cohort_log.txt"

* Settings. The defaults are sensible.
*
*   maxgroups     the most groups we will list for one column. Above this we
*                 just say how many there are, so we never end up printing a
*                 line for nearly every row in the data.
*   numgroups     a number column with more different values than this is
*                 treated as a measurement (like age) rather than as groups.
*   showextremes  show the smallest and largest value of each measurement.
*   dateparse     if this share of a text column reads as a date, call it a date.
*   iduniq        a whole-number column where nearly every row differs is
*                 almost certainly an identifier, so we hide what it contains.
*   hidesmall     set above 1 to hide groups with fewer than that many rows.
*   roundto       set above 1 to round the counts before printing them.
if "$maxgroups"    == "" global maxgroups    200
if "$numgroups"    == "" global numgroups    25
if "$showextremes" == "" global showextremes 1
if "$dateparse"    == "" global dateparse    0.90
if "$iduniq"       == "" global iduniq       0.98
if "$hidesmall"    == "" global hidesmall    5
if "$roundto"      == "" global roundto      1

* Save everything to a text file as well as showing it on screen.
capture log close _all
log using "$logfile", replace text name(profile)

display as text "{hline 78}"
display as text "Looking at: " as result `"`infile'"'
display as text "Date: " as result "$S_DATE $S_TIME"
if ${hidesmall} > 1 | ${roundto} > 1 {
    display as text "Small groups hidden, counts rounded."
}
else {
    display as text "Exact counts. Identifiers and long notes are not printed."
}
display as text "{hline 78}"

* Open the data and see how big it is.
use `"`infile'"', clear

local N = _N                    // rows, usually one per patient
quietly ds                      // ask Stata for the column names
local varlist `r(varlist)'
local K : word count `varlist'

display _newline as text "Rows: " as result `N' as text "   Columns: " as result `K'

if `N' == 0 {
    display as error "This file has no rows in it."
    log close profile
    exit
}


*==============================================================================
* A helper that works out what sort of column we are looking at.
*
* It gives back one of five answers:
*   identifier    a patient number, code or key. We never print the contents.
*   date          a date, however it happens to be stored.
*   measurement   something you would average, like age or a dose.
*   groups        a set of categories, like sex or stage or treating trust.
*   notes         long free text. We never print the contents.
*==============================================================================
capture program drop _classify
program define _classify, rclass
    args v

    quietly count
    local N = r(N)

    local stype : type `v'      // how Stata stores it, e.g. double or str10
    local vfmt  : format `v'    // how Stata shows it; dates begin with %t

    * Is this text rather than numbers?
    local istext = 0
    if substr("`stype'",1,3) == "str" local istext = 1

    * How many rows are blank?
    quietly count if missing(`v')
    local nblank = r(N)
    local nfilled = `N' - `nblank'

    * How many different values are there? We mark one row for each different
    * value, then count the marks. The second attempt is a fallback for unusual
    * column types where the first is not allowed.
    tempvar tag
    capture quietly egen `tag' = tag(`v') if !missing(`v')
    if _rc == 0 {
        quietly count if `tag' == 1
        local ndiff = r(N)
    }
    else {
        capture quietly bysort `v': gen byte `tag' = (_n == 1) if !missing(`v')
        if _rc == 0 {
            quietly count if `tag' == 1
            local ndiff = r(N)
        }
        else local ndiff = .
    }
    if missing(`ndiff') local ndiff = 0

    * Share of rows holding a different value. Near 1 means nearly all unique.
    local uniqshare = cond(`N' > 0, `ndiff'/`N', 0)

    * Does the name suggest it identifies a person or an address?
    local lname = lower("`v'")
    local nameid = 0
    foreach pat in patientid pseudo tumourid nhsnumber lsoa postcode _id {
        if strpos("`lname'","`pat'") local nameid = 1
    }
    if "`lname'" == "id" local nameid = 1

    * Does the name suggest an organisation code, such as a trust or hospital?
    * These have plenty of different values but they are still groups, and we
    * do want to see the full list of them.
    local nameorg = 0
    foreach pat in trust hosp provider site org procode alliance ccg {
        if strpos("`lname'","`pat'") local nameorg = 1
    }

    * Does the name suggest a clinical classification code, such as a
    * morphology or ICD code? These are whole numbers with a fair few different
    * values, but the number itself has no size or order to it - averaging a
    * morphology code makes no more sense than averaging a phone number. Treat
    * these the same way as organisation codes: always a set of groups.
    local namecode = 0
    foreach pat in morph icd_ icdcode _icd opcs {
        if strpos("`lname'","`pat'") local namecode = 1
    }
    local nameorg = max(`nameorg', `namecode')

    * Does it look like a date? Three clues, taken in turn: Stata already shows
    * it as a date, the name mentions a date, or the text actually reads as one.
    local isdate    = 0
    local textdate  = 0
    local dorder    = ""
    if substr("`vfmt'",1,2) == "%t"  local isdate = 1
    if strpos("`lname'","date")      local isdate = 1
    if strpos("`lname'","mdy")       local isdate = 1

    * A text column can hold dates whatever it is called, for instance
    * treat_start. So we try reading it three ways round and keep whichever
    * works for most rows. If nearly all rows read as a date, it is a date.
    if `istext' == 1 & `nfilled' > 0 {
        tempvar dayfirst monthfirst yearfirst
        capture quietly gen double `dayfirst'   = date(`v', "DMY")
        capture quietly gen double `monthfirst' = date(`v', "MDY")
        capture quietly gen double `yearfirst'  = date(`v', "YMD")
        if _rc == 0 {
            quietly count if !missing(`dayfirst')
            local okday = r(N)
            quietly count if !missing(`monthfirst')
            local okmonth = r(N)
            quietly count if !missing(`yearfirst')
            local okyear = r(N)

            local best = max(`okday', `okmonth', `okyear')
            if `best'/`nfilled' >= ${dateparse} {
                local isdate   = 1
                local textdate = 1
                if      `best' == `okday'   local dorder "DMY"
                else if `best' == `okmonth' local dorder "MDY"
                else                        local dorder "YMD"
            }
        }
    }

    * A number column made only of whole numbers, where nearly every row is
    * different, is almost certainly an identifier. Checking this stops a
    * patient number with an innocent name being summarised like an age.
    local wholenumbers = 0
    if `istext' == 0 & `nfilled' > 0 {
        capture quietly count if !missing(`v') & `v' != int(`v')
        if _rc == 0 & r(N) == 0 local wholenumbers = 1
    }

    * Long text is treated as notes and never printed. We judge this on how
    * varied the entries are and how long they are, rather than on a plain count
    * of different values, so that a trust or hospital code column is not
    * mistaken for notes.
    local typicallen = 0
    if `istext' == 1 & `nfilled' > 0 {
        tempvar len
        capture quietly gen `len' = length(`v') if !missing(`v')
        if _rc == 0 {
            quietly summarize `len', detail
            local typicallen = r(p50)
        }
    }

    * Text that is neither an organisation code nor a short tidy category, and
    * comes in many different varieties, is almost certainly notes.
    local looksnotes = 0
    if `istext' == 1 & `nameorg' == 0 {
        if `typicallen' > 40 & `ndiff' > 50                  local looksnotes = 1
        if `uniqshare' > 0.25 & `ndiff' > 100                local looksnotes = 1
    }

    * Now choose. The first rule that fits wins.
    *
    * Something is only an identifier if it really does hold a lot of different
    * values. Without that check a column such as
    * NHSE_reversed_imd_quintile_lsoas is hidden just because its name ends in
    * lsoas, and in a small file a column with a dozen categories can look
    * unique by accident. The last line catches a genuine id in a small file:
    * named like an id, and every row different.
    local manyvalues = (`ndiff' > ${numgroups})

    local class "groups"
    if      `isdate' == 1                                                         local class "date"
    else if `nameid' == 1 & `manyvalues' == 1                                     local class "identifier"
    else if `nameid' == 1 & `uniqshare' >= 0.99                                   local class "identifier"
    else if `istext' == 1 & `nameorg' == 0 & `manyvalues' == 1 & `uniqshare' > 0.5 local class "identifier"
    else if `istext' == 0 & `wholenumbers' == 1 & `manyvalues' == 1 & `uniqshare' > ${iduniq} local class "identifier"
    else if `looksnotes' == 1                                                     local class "notes"
    else if `istext' == 0 & `nameorg' == 0 & `manyvalues' == 1                    local class "measurement"

    * Organisation codes are always listed out in full, however many there are.
    local alwayslist = `nameorg'

    return local class    "`class'"
    return local dorder   "`dorder'"
    return scalar ndiff     = `ndiff'
    return scalar nblank    = `nblank'
    return scalar istext    = `istext'
    return scalar textdate  = `textdate'
    return scalar alwayslist = `alwayslist'
end


*==============================================================================
* Print a summary for one column, chosen to suit what it holds.
*==============================================================================
capture program drop _profilevar
program define _profilevar
    args v

    quietly count
    local N = r(N)

    * Ask the helper above what sort of column this is.
    _classify `v'
    local class    "`r(class)'"
    local dorder   "`r(dorder)'"
    local ndiff    = r(ndiff)
    local nblank   = r(nblank)
    local istext   = r(istext)
    local textdate = r(textdate)
    local alwayslist = r(alwayslist)

    local stype : type `v'
    local vfmt  : format `v'
    local vlab  : variable label `v'      // the plain description, if there is one
    local pblank = cond(`N' > 0, 100*`nblank'/`N', 0)
    local nfilled = `N' - `nblank'

    * Heading: the name, how it is stored, and what we think it holds.
    display _newline as result "`v'" as text "  [`stype', `class']"
    if `"`vlab'"' != "" display as text "    description: " as result `"`vlab'"'
    display as text "    blank: " as result %9.0f `nblank' ///
        as text " (" as result %4.1f `pblank' as text "%)   different values: " ///
        as result %6.0f `ndiff'

    * Nothing to say if every row is blank.
    if `nfilled' == 0 {
        display as text "    (every row is blank, so nothing to summarise)"
        exit
    }

    * Identifiers and notes: we say how long they are, never what they say.
    if inlist("`class'", "identifier", "notes") {
        if `istext' == 1 {
            tempvar len
            capture quietly gen `len' = length(`v') if !missing(`v')
            if _rc == 0 {
                quietly summarize `len', detail
                display as text "    contents not printed. typical length " ///
                    as result r(p50) as text " characters, longest around " ///
                    as result r(p95)
            }
            else display as text "    contents not printed."
        }
        else {
            display as text "    a number used as an identifier. contents not printed."
        }
        exit
    }

    * Dates: show the span and the middle date, not the individual dates.
    if "`class'" == "date" {
        if `textdate' == 1 {
            * Stored as text, so make a proper date copy before summarising.
            tempvar realdate
            quietly gen double `realdate' = date(`v', "`dorder'")
            format `realdate' %td
            quietly summarize `realdate', detail
            display as text "    stored as text, read as `dorder' order."
            display as text "    runs from " as result %td r(p1) ///
                as text " to " as result %td r(p99) ///
                as text ", middle date " as result %td r(p50)

            * Flag any rows that could not be read as a date at all.
            quietly count if missing(`realdate') & !missing(`v')
            if r(N) > 0 display as text "    note: " as result r(N) ///
                as text " rows could not be read as a date"
        }
        else if `istext' == 1 {
            display as text "    named like a date, but the contents are not dates." ///
                " not printed."
        }
        else {
            quietly summarize `v', detail
            if substr("`vfmt'",1,2) == "%t" {
                display as text "    runs from " as result %td r(p1) ///
                    as text " to " as result %td r(p99) ///
                    as text ", middle date " as result %td r(p50)
            }
            else {
                display as text "    named like a date but held as a plain number." ///
                    " middle value " as result r(p50)
            }
        }
        exit
    }

    * Measurements: the average, the spread, and where most rows sit.
    if "`class'" == "measurement" {
        quietly summarize `v', detail
        display as text "    average " as result %10.2f r(mean) ///
            as text "   spread " as result %10.2f r(sd)
        if ${showextremes} == 1 {
            display as text "    smallest " as result r(min) ///
                as text "   largest " as result r(max)
        }

        * The values below which 1, 5, 25, 50, 75, 95 and 99 percent of rows sit.
        display as text "    1% / 5% / 25% / half / 75% / 95% / 99% of rows fall below:"
        display as text "      " ///
            as result r(p1)  as text " / " as result r(p5)  as text " / " ///
            as result r(p25) as text " / " as result r(p50) as text " / " ///
            as result r(p75) as text " / " as result r(p95) as text " / " ///
            as result r(p99)

        * If the average sits well above the middle value, a few big values are
        * dragging it up. Worth knowing before anyone quotes the average.
        local lopsided = r(mean) - r(p50)
        display as text "    (average minus middle value = " as result %6.2f `lopsided' ///
            as text ". a positive figure means a few large values pull the average up)"

        * Negatives often mean a data problem, such as two dates the wrong way round.
        quietly count if `v' < 0 & !missing(`v')
        if r(N) > 0 display as text "    note: " as result r(N) ///
            as text " rows hold a negative value, worth a look"
        exit
    }

    * Everything else is a set of groups, so we tabulate it.

    * Numbers often come with a lookup giving each code a meaning, for example
    * 1 means Male. See whether this column has one.
    local lookup : value label `v'
    if "`lookup'" != "" {
        display as text "    uses the lookup called " as result "`lookup'"
    }
    else if `istext' == 0 {
        display as text "    (no lookup attached, so the raw numbers are shown)"
    }

    * Some codes have no entry in the lookup at all. To spot those we make a
    * translated copy of the column: any row that fails to translate is a code
    * with no meaning attached. We cannot simply compare the meaning against
    * the number, because a code of 3 is often labelled "3".
    tempvar translated
    local havelookup = 0
    if "`lookup'" != "" & `istext' == 0 {
        capture decode `v', generate(`translated')
        if _rc == 0 local havelookup = 1
    }

    * Get the list of groups. There is deliberately no clean option here: it
    * would strip the quote marks, and a group called "Emergency presentation"
    * would then be split at the space into two groups that do not exist.
    quietly levelsof `v', local(levels)
    local ngroups : word count `levels'

    * If there are so many groups that printing them would amount to listing
    * the column row by row, we stop and just report how many there are.
    * Trust and hospital columns are always listed out, however many there are,
    * because seeing the full set of codes is usually the point.
    if `ngroups' > ${maxgroups} & `alwayslist' == 0 {
        display as text "    " as result `ngroups' ///
            as text " different groups, too many to list here."
        exit
    }

    * First pass: count the rows in each group and remember the counts.
    local i = 0
    foreach lev of local levels {
        local ++i
        if `istext' == 1 quietly count if `v' == `"`lev'"'
        else             quietly count if `v' == `lev'
        local cnt`i' = r(N)
        local lev`i' `"`lev'"'
    }

    * Decide whether to hide any small groups. Off unless you asked for it.
    * Hiding exactly one group would let someone work out its size by taking
    * the others away from the total, so in that case we hide a second one.
    local hiding   = (${hidesmall} > 1)
    local rounding = (${roundto} > 1)

    forvalues j = 1/`ngroups' {
        local hide`j' = 0
    }
    local hidsecond = 0
    if `hiding' {
        local nhidden = 0
        forvalues j = 1/`ngroups' {
            local hide`j' = (`cnt`j'' < ${hidesmall})
            if `hide`j'' local ++nhidden
        }
        if `nhidden' == 1 & `ngroups' > 1 {
            * find the smallest group still showing and hide that one too
            local smallest = 0
            local smallestn = .
            forvalues j = 1/`ngroups' {
                if !`hide`j'' & `cnt`j'' < `smallestn' {
                    local smallestn = `cnt`j''
                    local smallest = `j'
                }
            }
            if `smallest' > 0 {
                local hide`smallest' = 1
                local hidsecond = 1
            }
        }
    }

    * The figure the percentages are worked out of.
    if `rounding' local base = round(`nfilled', ${roundto})
    else          local base = `nfilled'

    if `hiding' | `rounding' {
        display as text "    groups (out of " as result `base' ///
            as text " rows with a value; small groups hidden, counts rounded):"
        if `hidsecond' == 1 {
            display as text "      [a second group is hidden so the first cannot" ///
                " be worked out by subtraction]"
        }
    }
    else {
        display as text "    groups (out of " as result `base' as text " rows with a value):"
    }

    * Second pass: print each group, with its meaning where we have one.
    forvalues j = 1/`ngroups' {
        local thislev `"`lev`j''"'
        local thiscnt = `cnt`j''

        if `istext' == 1 {
            local shown `"`thislev'"'
        }
        else {
            local shown "`thislev'"
            if `havelookup' == 1 {
                * rows that failed to translate have no entry in the lookup
                quietly count if `v' == `thislev' & missing(`translated')
                if r(N) > 0 {
                    local shown "`thislev' = (no meaning attached to this code)"
                }
                else {
                    local meaning : label `lookup' `thislev'
                    local shown "`thislev' = `meaning'"
                }
            }
            else if "`lookup'" != "" {
                local meaning : label `lookup' `thislev'
                local shown "`thislev' = `meaning'"
            }
        }

        if `hide`j'' == 1 {
            local printn "hidden"
            local printpc "."
        }
        else {
            if `rounding' local shownn = round(`thiscnt', ${roundto})
            else          local shownn = `thiscnt'
            local printn "`shownn'"
            local printpc = string(100*`shownn'/`base', "%4.1f")
        }
        display as text "      " as result `"`shown'"' ///
            as text "   n=" as result "`printn'" ///
            as text "   (" as result "`printpc'" as text "%)"
    }
end


* Go through every column. If one causes trouble we note it and carry on,
* rather than stopping the whole run.
display _newline as text "{hline 78}"
display as text "1. Each column in turn"
display as text "{hline 78}"

local failed ""
foreach v of local varlist {
    capture noisily _profilevar `v'
    if _rc != 0 {
        display as error "    [could not summarise `v', carrying on]"
        local failed "`failed' `v'"
    }
}
if "`failed'" != "" {
    display _newline as error "Columns that could not be summarised:`failed'"
}


* The lookups in full.
*
* The tables above show only the codes that actually turn up. A lookup can also
* define codes that never appear, such as "9 = Not known" in a file where
* nobody is coded 9. Those still tell you how the data was meant to be coded,
* so every lookup is listed here in full.
display _newline as text "{hline 78}"
display as text "2. The lookups, and what each code means"
display as text "{hline 78}"

quietly label dir                  // the names of all the lookups in the file
local alllookups `r(names)'

if `"`alllookups'"' == "" {
    display as text "  (this file has no lookups, so codes are shown as they are)"
}
else {
    foreach L of local alllookups {
        * which columns use this lookup?
        local users ""
        foreach v of local varlist {
            local thislookup : value label `v'
            if "`thislookup'" == "`L'" local users "`users' `v'"
        }
        if "`users'" == "" local users " (not used by any column)"
        display _newline as result "  `L'" as text "   used by:" as result "`users'"
        capture noisily label list `L'
    }
}


* A single line per column, handy to copy into notes or an email.
display _newline as text "{hline 92}"
display as text "3. Summary table (column | held as | what it is | different values | % blank | lookup)"
display as text "{hline 92}"

foreach v of local varlist {
    capture _classify `v'
    if _rc != 0 {
        display as text "  " as result %-30s "`v'" as text " | " ///
            as error "could not work this one out"
        continue
    }
    local class  "`r(class)'"
    local ndiff  = r(ndiff)
    local nblank = r(nblank)
    local pblank = cond(`N' > 0, 100*`nblank'/`N', 0)

    local stype  : type `v'
    local lookup : value label `v'
    local lookupshown = cond("`lookup'" == "", "-", "`lookup'")

    display as text "  " as result %-30s "`v'" as text " | " ///
        as result %-9s "`stype'" as text " | " as result %-11s "`class'" ///
        as text " | " as result %6.0f `ndiff' as text " | " ///
        as result %5.1f `pblank' as text "% | " as result "`lookupshown'"
}

display _newline as text "{hline 78}"
display as text "Finished. The same write-up has been saved to:"
display as text "  $logfile"
display as text "{hline 78}"

log close profile
