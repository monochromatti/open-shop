// Article layout, using only fonts bundled with Typst.
#let ink = rgb("#20262a")
#let accent = rgb("#214c62")
#let muted = rgb("#56636b")
#let rule = rgb("#bac5cb")

#let article(body) = {
  set document(title: "Nonlinear hydropower scheduling with coupled hydraulics and delayed river networks")
  set page(
    paper: "a4",
    margin: (left: 36mm, right: 22mm, top: 23mm, bottom: 23mm),
    header: context {
      if counter(page).get().first() > 1 {
        set text(size: 8.5pt, fill: muted)
        block(width: 100%, inset: (bottom: 5pt), stroke: (bottom: 0.4pt + rule))[
          Nonlinear hydropower scheduling
        ]
      }
    },
    footer: context align(right,
      text(size: 8.5pt, fill: muted)[
        #counter(page).display() / #counter(page).final().first()
      ]),
  )
  set text(font: "Libertinus Serif", size: 10.5pt, fill: ink, lang: "en")
  set par(justify: true, leading: 0.6em, spacing: 0.65em, first-line-indent: 1em)
  set heading(numbering: "1.1")
  set math.equation(numbering: "(1)", supplement: none)
  show math.equation: set text(font: "New Computer Modern Math")
  show math.equation.where(block: true): set block(above: 0.85em, below: 0.85em)
  show heading.where(level: 1): it => {
    set text(size: 13.5pt, weight: "bold", fill: accent)
    set block(above: 1.6em, below: 0.7em)
    it
  }
  show heading.where(level: 2): it => {
    set text(size: 11pt, weight: "bold")
    set block(above: 1.15em, below: 0.55em)
    it
  }
  show figure.caption: it => {
    set text(size: 9pt, fill: muted)
    set par(justify: false, leading: 0.5em, first-line-indent: 0pt)
    align(left, it)
  }
  show figure.where(kind: table): set figure.caption(position: top)
  show table: set text(size: 9.5pt)
  show link: set text(fill: accent)
  show ref: it => {
    set text(fill: accent)
    if it.element != none and it.element.func() == math.equation {
      link(it.target, numbering("(1)", ..counter(math.equation).at(it.element.location())))
    } else {
      it
    }
  }
  show bibliography: it => {
    set text(size: 9pt)
    set par(leading: 0.5em, spacing: 0.5em)
    it
  }
  body
}

#let title-block(body) = block(below: 20pt)[
  #set par(justify: false, leading: 0.5em, first-line-indent: 0pt)
  #text(size: 26pt, fill: accent, body)
]

#let guide-entry(target, title) = block(above: 10pt)[
  #link(target, text(size: 9pt, weight: "bold", title)) \
  #context text(size: 8.5pt, fill: muted)[Page #counter(page).at(target).first()]
]

#let abstract-block(body) = grid(
  columns: (29mm, 1fr),
  gutter: 7mm,
  {
    set par(justify: false, leading: 0.5em, first-line-indent: 0pt)
    text(size: 8.5pt, fill: muted, smallcaps[In this article])
    guide-entry(<sec:water>, [Physical model])
    guide-entry(<sec:solution>, [Solution and checks])
    guide-entry(<sec:graphs>, [Curve equations and bounds])
  },
  block(inset: (top: 9pt, bottom: 10pt),
    stroke: (top: 0.6pt + rule, bottom: 0.6pt + rule))[
    #set text(size: 10pt)
    #set par(leading: 0.55em, first-line-indent: 0pt)
    #text(size: 9pt, weight: "bold", fill: accent, smallcaps[Abstract])
    #parbreak()
    #body
  ],
)
