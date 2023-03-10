(cl:in-package #:common-lisp-user)

(defpackage #:serial
  (:use #:cl)
  (:export #:serializable))


(defpackage #:topology
  (:use #:common-lisp)
  (:nicknames #:ts)
  (:export
   #:make-topology-from-residue
   #:topologyp
   #:degree-difference
   #:constitution-atom-index
   #:atom-name
   #:element
   #:bonds
   #:index
   #:to-atom-index
   #:order
   #:children
   #:name
   #:plug-bonds
   #:atom-name
   #:bond-order
   #:plug-bond
   #:root-node
   #:constitution-from-graph
   #:stereo-information
   #:stereoisomer
   #:make-in-plug
   #:make-out-plug
   #:add-plug
   #:plugs
   #:name
   #:configurations
   #:make-constitution-atoms-from-residue
   #:topology-list
   #:property-list
   #:topology
   #:in-plug
   #:has-plug-named
   #:plug-named
   #:has-in-coupling-p
   #:define-topology
   #:constitution
   #:constitution-atoms
   #:stereoisomer-name
   #:stereoisomer-atoms
   #:build-residue-single-name
   #:connect-residues
   #:is-in-plug-name
   #:other-plug-name
   #:find-in-plug
   #:atom-type
   #:joint-template
   #:bonded-joint-template
   #:in-plug-bonded-joint-template
   #:complex-bonded-joint-template
   #:jump-joint-template
   #:in-plug
   #:name
   #:id
   #:source-monomer
   #:target-monomer
   #:source-plug-name
   #:target-plug-name
   #:children
   #:constitution
   #:atom-name
   #:parent
   #:input-stub-joints
   #:plugs
   #:add-child
   #:sibling
   #:aggregate
   #:ordered-monomers

   #:seen-fragment-internals
   #:bad-fragment-internals
   #:load-fragment-conformations
   #:fragment-conformations
   #:fragment-internals
   #:out-of-focus-internal
   #:out-of-focus-internals
   #:internal
   #:jump-internal
   #:bonded-internal
   #:complex-bonded-internal
   #:fragments
   #:internals
   #:bond
   #:angle
   #:dihedral
   #:foldamer-monomer-context
   #:save-fragment-conformations
   #:load-fragment-conformations
   #:fragment-conformations
   #:fragment-conformations-map
   #:monomer-context-to-fragment-conformations
   #:dump-fragment-internals
   #:index
   #:total-count
   #:bad-geometry-p
   
   #:make-in-plug-bonded-joint-template
   #:make-bonded-joint-template
   #:make-jump-joint-template
   #:make-complex-bonded-joint-template
   #:stereo-information
   #:build-internal-coordinate-joint-template-tree
   #:build-joint-template
   #:extract-prepare-topologys
   #:target-monomer
   #:define-cap
   #:build-molecule
   #:oligomer-space
   #:oligomers
   #:monomer
   #:directional-coupling
   #:out-plug-name
   #:in-plug-name
   #:monomers
   #:couplings
   #:number-of-sequences
   #:current-stereoisomer-name
   #:other-monomer
   #:goto-sequence
   #:monomer-positions
   #:ataggregate
   #:atmolecules
   #:atresidues
   #:joints
   #:zero-all-atom-tree-external-coordinates
   #:build-all-atom-tree-external-coordinates
   #:oligomer-space-directional-coupling-iterator-factory

   #:make-joint-tree
   #:joint-tree
   #:children
   #:root
   #:atom-id-map
   #:make-conformation
   #:make-oligomer
   #:topologys-in-oligomer-space

   #:in-plug-name-p
   #:out-plug-name-p
   #:coupling-name
   
   ;;;
   #:update-joint-tree-internal-coordinates
   #:build-all-atom-tree-external-coordinates
   #:define-part
   #:*topology-groups*
   #:restraints

   #:missing-fragment-match
   #:matched-fragment-conformations-map
   #:fragment-matches
   #:oligomer-monomer-name-at-index
   #:oligomer-monomer-name-for-monomer
   #:missing-fragment-matches-count
   #:matched-fragment-conformations-summary

   #:fragment-match-key-before-monomer-context-index
   #:fragment-match-key-after-monomer-context-index
   #:missing-fragment-match-key-before-monomer-context-index
   #:missing-fragment-match-key-after-monomer-context-index
   #:monomer-contexts-vector
   #:focus-monomer-name
   
   #:build-shape
   #:monomer-shapes
   #:monomer-shape-vector
   #:root-monomer
   #:random-fragment-conformation-index
   #:make-oligomer-shape
   #:build-one-molecule-for-topology
   #:sketch-svg
   #:make-oligomer-space
   #:calculate-number-of-sequences
   #:energy-function
   #:copy-fragment-internals
   #:copy-internal
   #:deg-to-rad
   #:rad-to-deg
   #:radians-add
   #:radians-incf
   #:radians-difference
   #:angle-sub
   #:walk-joint-template
   #:properties
   #:degrees-add
   #:probability
   #:oligomer-force-field-name
   #:oligomer
   #:build-residue-for-isomer
   #:stereoisomer-atom
   #:atom-charge))

(defpackage #:monomer-context
  (:use #:common-lisp)
  (:export
   #:parse
   #:matches
   #:match-as-symbol
   #:match
   #:match-iterator
   ))

(defpackage #:topology.graphviz
  (:use #:common-lisp)
  (:nicknames #:topg)
  (:export
   #:nodes
   #:node-label
   #:directed-edges
   #:directed-edge-from
   #:directed-edge-to
   #:directed-edge-label
   #:render-dag
   #:undirected-edges
   #:undirected-edge-from
   #:undirected-edge-to
   #:undirected-edge-label
   #:graph-label
   #:graph-name
   #:render-foldamer-joint-trees
   #:make-graph
   #:node-id
   #:dot-svg-foldamer-joint-trees))


