(in-package :topology)


(defclass monomer-shape ()
  ((fragment-conformation-index :initarg :fragment-conformation-index)
   (monomer :initarg :monomer :accessor monomer)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (monomer-context-index :initarg :monomer-context-index :accessor monomer-context-index)))

(defclass oligomer-shape ()
  ((oligomer :initarg :oligomer :accessor oligomer)
   (matched-fragment-conformations-map :initarg :matched-fragment-conformations-map :accessor matched-fragment-conformations-map)
   (foldamer :initarg :foldamer :accessor foldamer)
   (monomer-shape-vector :initarg :monomer-shape-vector :accessor monomer-shape-vector)
   (monomer-shape-map :initarg :monomer-shape-map :accessor monomer-shape-map)
   (the-root-monomer :initarg :the-root-monomer :accessor the-root-monomer)
   (in-monomers :initarg :in-monomers :accessor in-monomers)
   (out-monomers :initarg :out-monomers :accessor out-monomers)
   ))


(defun make-oligomer-shape (oligomer matched-fragment-conformations-map foldamer)
  (multiple-value-bind (monomer-shape-vector the-root-monomer in-monomers out-monomers monomer-shape-map)
      (loop with monomer-shape-vector = (make-array (length (monomers oligomer)))
            with in-monomers = (make-hash-table)
            with out-monomers = (make-hash-table)
            with the-root-monomer = nil
            with monomer-shape-map = (make-hash-table)
            for index from 0
            for monomer across (monomers oligomer)
            for monomer-context = (topology:foldamer-monomer-context monomer oligomer foldamer)
            for monomer-context-index = (gethash monomer-context (monomer-context-index-map matched-fragment-conformations-map))
            for monomer-shape = (make-instance 'monomer-shape
                                               :monomer monomer
                                               :monomer-context monomer-context
                                               :monomer-context-index monomer-context-index)
            for couplings = (couplings monomer)
            for in-monomer = nil
            for out-mons = nil
            do (format t "monomer = ~a~%" monomer)
            do (setf (gethash monomer monomer-shape-map) monomer-shape)
            do (maphash (lambda (key coupling)
                          (if (in-plug-name-p key)
                              (progn
                                (setf in-monomer (topology:source-monomer coupling))
                                (setf (gethash monomer in-monomers) (topology:source-monomer coupling))
                                (format t "In plug coupling ~a ~a~%" key coupling))
                              (progn
                                (push (topology:target-monomer coupling) out-mons)
                                (format t "Out plug coupling ~a ~a~%" key coupling))))
                        couplings)
            do (unless in-monomer
                 (setf the-root-monomer monomer))
            do (setf (gethash monomer out-monomers) out-mons)
            do (setf (aref monomer-shape-vector index) monomer-shape)
            do (format t "monomer-context ~a~%" monomer-context)
            finally (return (values monomer-shape-vector the-root-monomer in-monomers out-monomers monomer-shape-map)))
    (make-instance 'oligomer-shape
                   :oligomer oligomer
                   :matched-fragment-conformations-map matched-fragment-conformations-map
                   :foldamer foldamer
                   :monomer-shape-vector monomer-shape-vector
                   :monomer-shape-map monomer-shape-map
                   :the-root-monomer the-root-monomer
                   :in-monomers in-monomers
                   :out-monomers out-monomers)))


(defun all-monomers-impl (root shape)
  (format t "monomer ~a in: ~a~%" root (gethash root (in-monomers shape)))
  (let ((out-monomers (gethash root (out-monomers shape))))
    (loop for out-monomer in out-monomers
          do (all-monomers-repl out-monomer shape))))

(defun all-monomers (shape)
  (let ((root (the-root-monomer shape)))
    (all-monomers-impl root shape)))




(defun random-fragment-conformation-index-impl (root shape)
  (format t "monomer ~a in: ~a~%" root (gethash root (in-monomers shape)))
  (let ((out-monomers (gethash root (out-monomers shape))))
    (loop for out-monomer in out-monomers
          do (all-monomers-repl out-monomer shape))))

(defun random-fragment-conformation-index (shape)
  (let ((root (the-root-monomer shape)))
    (random-fragment-conformation-index-impl root shape)))
