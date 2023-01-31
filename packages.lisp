(cl:in-package #:common-lisp-user)

(defpackage #:serial
  (:use #:cl)
  (:export #:serializable))


(defpackage #:topology
  (:use #:common-lisp)
  (:nicknames #:ts)
  (:export
   #:constitution-atom-index
   #:atom-name
   #:element
   #:bonds
   #:index
   #:to-atom-index
   #:order
   #:children
   #:name
   #:atom-names
   #:root-node
   #:constitution-from-graph
   #:stereo-information
   #:plugs
   #:name
   #:configurations
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
   #:oligomer
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
   #:copy-joint-positions-into-atoms
   #:copy-atom-positions-into-joints
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
   #:sketch-svg))

(defpackage #:monomer-context
  (:use #:common-lisp)
  (:export
   #:parse
   #:matches
   #:match-as-symbol
   #:match
   ))

(defpackage #:foldamer
  (:use #:common-lisp)
  (:nicknames #:fd)
  (:export
   #:define-foldamer
   #:parse
   #:oligomer-space
   #:unused-trainer-contexts
   #:maybe-remove-unused-trainers
   #:generate-training-oligomers
   #:find-oligomer-for-monomer-context
   #:oligomer-monomer-context-iterator
   #:build-trainer
   #:prepare-to-build-trainer
   #:foldamer-monomer-context
   #:calculate-files
   #:verify-foldamer-describes-oligomer-space
   #:valid-trainer-contexts
   #:extract-fragment-conformations-map

   #:foldamer-setup
   #:foldamer-status
   #:foldamer-extract-conformations
   #:foldamer-run-node

   #:trainer-context
   #:trainer-job
   #:node-index
   #:optimize-fragment-conformations-map
   #:foldamer-describe-missing-fragment-matches
   #:monomer-contents-vector
   #:foldamer-describe-missing-match
   #:angle-difference
   #:monomer-context
   #:monomer-context-to-oligomer-map
   #:topologys

   #:verify-all-training-molecules-can-be-parameterized
   #:load-force-field
   #:load-foldamer-conformations-map
   #:save-foldamer-conformations-map))


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
