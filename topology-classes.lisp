(in-package :topology)

(defclass constitution-atom ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (index :initarg :index :accessor index)
   (element :initarg :element :accessor element)
   (atom-type :initarg :atom-type :accessor atom-type)
   (stereochemistry-type :initarg :stereochemistry-type :accessor stereochemistry-type)
   (properties :initform nil :initarg :properties :accessor properties)
   (bonds :initform nil :initarg :bonds :accessor bonds)))

(cando:make-class-save-load constitution-atom
                            :print-unreadably
                            (lambda (obj stream)
                              (print-unreadable-object (obj stream :type t)
                                (format stream "~a" (atom-name obj)))))

(defclass constitution-bond ()
  ((to-atom-index :initarg :to-atom-index :accessor to-atom-index)
   (order :initarg :order :accessor order)))


(cando:make-class-save-load constitution-bond)

(defclass constitution-atoms ()
  ((atoms :initarg :atoms :accessor atoms)))

(cando:make-class-save-load constitution-atoms)

(defclass stereoconfiguration ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (configuration :initarg :configuration :accessor configuration)))

(cando:make-class-save-load stereoconfiguration
                            :print-unreadably
                            (lambda (obj stream)
                              (print-unreadable-object (obj stream :type t)
                                (format stream "~a ~a" (atom-name obj) (configuration obj)))))

(defclass stereoisomer ()
  ((name :initarg :name :accessor name)
   (pdb :initarg :pdb :accessor pdb)
   (stereoisomer-index :initarg :stereoisomer-index :accessor stereoisomer-index)
   (configurations :initarg :configurations :accessor configurations)))

(cando:make-class-save-load stereoisomer
                            :print-unreadably
                            (lambda (obj stream)
                              (print-unreadable-object (obj stream :type t)
                                (format stream "~a ~s" (name obj) (configurations obj)))))

(defclass constitution ()
  ((name :initarg :name :accessor name)
   (constitution-atoms :initarg :constitution-atoms :accessor constitution-atoms)
   (plugs :type hash-table :initarg :plugs :accessor plugs)
   (topology-list :initform nil :initarg :topology-list :accessor topology-list)
   (stereo-information :initarg :stereo-information :accessor stereo-information)))

(cando:make-class-save-load constitution)

(defclass stereoisomer-atom ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (constitution-atom-index :initarg :constitution-atom-index :accessor constitution-atom-index)
   (atom-charge :initarg :atom-charge :accessor atom-charge)
   (atom-type :initarg :atom-type :accessor atom-type)))

(cando:make-class-save-load stereoisomer-atom)

(defclass stereoisomer-virtual-atom (stereoisomer-atom)
  ())


(defclass restraint ()
  ())

(defclass dihedral-restraint ()
  ((atom1-name :initarg :atom1-name :accessor atom1-name)
   (atom2-name :initarg :atom2-name :accessor atom2-name)
   (atom3-name :initarg :atom3-name :accessor atom3-name)
   (atom4-name :initarg :atom4-name :accessor atom4-name)
   (dihedral-min-degrees :initarg :dihedral-min-degrees
                         :accessor dihedral-min-degrees)
   (dihedral-max-degrees :initarg :dihedral-max-degrees
                         :accessor dihedral-max-degrees)
   (weight :initarg :weight :accessor weight)))

(cando:make-class-save-load dihedral-restraint
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a ~a ~a ~a ~a ~a"
             (atom1-name obj)
             (atom2-name obj)
             (atom3-name obj)
             (atom4-name obj)
             (dihedral-min-degrees obj)
             (dihedral-max-degrees obj)
             (weight obj)))))

(defclass topology ()
  ((name :initarg :name :accessor name)
   (constitution :initarg :constitution :accessor constitution)
   (property-list :initform nil :initarg :property-list :accessor property-list)
   (plugs :type hash-table :initarg :plugs :accessor plugs)
   (joint-template :initarg :joint-template :accessor joint-template)
   (restraints :initform nil :initarg :restraints :accessor restraints)
   (stereoisomer-atoms :initform (make-array 4 :adjustable t :fill-pointer 0)
                       :initarg :stereoisomer-atoms :accessor stereoisomer-atoms)
   ))


(cando:make-class-save-load topology
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (name obj)))))

(defun has-plug-named (topology plug-name)
  (gethash plug-name (plugs topology)))

(defun plug-named (topology plug-name)
  (gethash plug-name (plugs topology)))

(defun find-in-plug (topology)
  (maphash (lambda (name plug)
             (declare (ignore name))
             (when (typep plug 'in-plug)
               (return-from find-in-plug plug)))
           (plugs topology)))

(defun out-plugs-as-list (topology)
  (let (out-plugs)
    (maphash (lambda (name plug)
               (declare (ignore name))
               (when (typep plug 'out-plug)
                 (push plug out-plugs)))
             (plugs topology))
    out-plugs))

(defun all-out-plug-names-that-match-in-plug-name (topology in-plug-name)
  (let (out-plugs
        (in-coupling-name (coupling-name in-plug-name)))
    (maphash (lambda (name plug)
               (declare (ignore name))
               (when (and (typep plug 'out-plug)
                          (eq (coupling-name (name plug))
                              in-coupling-name))
                 (push (name plug) out-plugs)))
             (plugs topology))
    out-plugs))

(defun plugs-as-list (topology)
  (let (plugs)
    (maphash (lambda (name plug)
               (declare (ignore name))
               (push plug plugs))
             (plugs topology))
    plugs))

(defclass plug ()
  ((name :initarg :name :accessor name)
   (atom-names :initform (make-array 16)
               :initarg :atom-names
               :accessor atom-names)))

(cando:make-class-save-load plug
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (name obj) (atom-names obj)))))

(defclass in-plug (plug)
  ())

(cando:make-class-save-load in-plug)

(defclass out-plug (plug)
  ())

(cando:make-class-save-load out-plug)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Monomer, coupling, and oligomer
;;;


(defclass monomer ()
  ((id :initform nil
       :initarg :id
       :accessor id)
   (couplings :type hash-table
              :initform (make-hash-table)
              :initarg :couplings
              :accessor couplings)
   (monomers :initarg :monomers :accessor monomers)))

(cando:make-class-save-load monomer
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream ":id ~a :monomers ~a" (id obj) (monomers obj)))))

(defun number-of-stereoisomers (monomer)
  (length (monomers monomer)))

(defun current-stereoisomer-name (monomer oligomer)
  (let ((monomer-index (position monomer (monomers (oligomer-space oligomer)))))
    (unless monomer-index
      (error "Could not find monomer ~a in ~a" monomer oligomer))
    (elt (monomers monomer) (aref (monomer-indices oligomer) monomer-index))))

(defun monomer-topologys-all-have-plug-named (monomer plug-name)
  (let ((topologys (loop for mon-name in (monomers monomer)
                         collect (chem:find-topology mon-name t))))
    (cond
      ((every (lambda (top) (has-plug-named top plug-name)) topologys)
       t)
      ((every (lambda (top) (not (has-plug-named top plug-name))) topologys)
       nil)
      (t (error "Some topologys of ~s have the plug named ~a and some do not" (monomers monomer) plug-name)))))

(defun monomer-topology (monomer oligomer)
  (let ((monomer-position (position monomer (monomers (oligomer-space oligomer)))))
    (unless monomer-position
      (error "Could not find ~a in ~a" monomer oligomer))
    (let* ((monomer (aref (monomers (oligomer-space oligomer)) monomer-position))
           (monomer-index (aref (monomer-indices oligomer) monomer-position))
           (monomer-name (elt (monomers monomer) monomer-index)))
      (chem:find-topology monomer-name t))))

(defun has-in-coupling-p (monomer)
  (maphash (lambda (plug-name coupling)
             (declare (ignore coupling))
             (when (in-plug-name-p plug-name)
               (return-from has-in-coupling-p t)))
           (couplings monomer)))

(defun monomer-plug-named (monomer plug-name)
  (gethash plug-name (couplings monomer)))

(defgeneric other-monomer (coupling monomer))

(defun monomer-on-other-side (monomer plug-name)
  (let ((coupling (monomer-plug-named monomer plug-name)))
    (unless coupling
      (error "Coupling for plug ~a was not found" plug-name))
    (other-monomer coupling monomer)))


(defclass coupling ()
  ((name :initarg :name :accessor name)))

(defclass directional-coupling (coupling)
  ((source-plug-name :initarg :source-plug-name :accessor source-plug-name)
   (target-plug-name :initarg :target-plug-name :accessor target-plug-name)
   (source-monomer :initarg :source-monomer :accessor source-monomer)
   (target-monomer :initarg :target-monomer :accessor target-monomer)))

(cando:make-class-save-load directional-coupling
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a ~a ~a" (source-monomer obj) (source-plug-name obj) (target-plug-name obj) (target-monomer obj)))))

(defmethod other-monomer ((coupling directional-coupling) monomer)
  (cond
    ((eq (source-monomer coupling) monomer)
     (target-monomer coupling))
    ((eq (target-monomer coupling) monomer)
     (source-monomer coupling))
    (t (error "Could not find monomer ~a in coupling ~a" monomer coupling))))


(defclass ring-coupling (coupling)
  ((plug1 :initarg :plug1 :accessor plug1)
   (plug2 :initarg :plug2 :accessor plug2)
   (monomer1 :initarg :monomer1 :accessor monomer1)
   (monomer2 :initarg :monomer2 :accessor monomer2)))

(defmethod other-monomer ((coupling ring-coupling) monomer)
  (cond
    ((eq (monomer1 coupling) monomer)
     (monomer2 coupling))
    ((eq (monomer2 coupling) monomer)
     (monomer1 coupling))
    (t (error "Could not find monomer ~a in coupling ~a" monomer coupling))))


(defun has-ring-closing-coupling (monomer)
  (maphash (lambda (plug-name coupling)
             (declare (ignore plug-name))
             (when (typep coupling 'ring-coupling)
               (return-from has-ring-closing-coupling t)))
           (couplings monomer)))

(defclass oligomer-space ()
  ((foldamer :initarg :foldamer :accessor foldamer)
   (monomers :initform (make-array 16 :adjustable t :fill-pointer 0)
             :initarg :monomers :accessor monomers)
   (couplings :initform (make-array 16 :adjustable t :fill-pointer 0)
              :initarg :couplings :accessor couplings)
   (%number-of-sequences :initform nil
                        :initarg :number-of-sequences
                        :accessor %number-of-sequences)))

(cando:make-class-save-load oligomer-space
  :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t))))


(defun topologys-in-oligomer-space (oligomer-space)
  (loop for monomer across (monomers oligomer-space)
        for monomer-names = (monomers monomer)
        append (loop for monomer-name in monomer-names
                     for topology = (chem:find-topology monomer-name t)
                     collect topology)))

(defun number-of-monomers (oligomer-space)
  (length (monomers oligomer-space)))

(defun add-monomer (oligomer-space monomer)
  (vector-push-extend monomer (monomers oligomer-space)))

(defun calculate-number-of-sequences (oligomer-space)
  (let ((num 1))
  (loop for monomer across (monomers oligomer-space)
        do (setf num (* num (length (monomers monomer)))))
    num))

(defun number-of-sequences (oligomer-space)
  ;; Lazy calculate number-of-sequences
  (when (null (%number-of-sequences oligomer-space))
    (setf (%number-of-sequences oligomer-space) (calculate-number-of-sequences oligomers-space)))
  (%number-of-sequences oligomer-space))

(defun make-oligomer-space (foldamer tree &key (parts *parts*))
  "Make an oligomer-space from a description in the **tree**.
The tree is a nested list of lists that look like 
(component coupling component coupling component ... ).
It starts with the name of a topology or a design:part.
The component can be a name of a single topology, it can be a design:part,
or it can 

Examples:
(make-oligomer-space '(:ccap :default :ala :default :ala :default :ser :default :ncap))
(make-oligomer-space '((:lego3 :label :first)
                 :default :lego3
                 :default :lego3
                 :default (cycle :first)))
"
  (let* ((oligomer-space (make-instance 'oligomer-space
                                        :foldamer foldamer))
         (labels (make-hash-table)))
    (interpret-subtree oligomer-space tree labels :parts parts)
    (setf (number-of-sequences oligomer-space)
          (calculate-number-of-sequences oligomer-space))
    oligomer-space))

(defclass oligomer ()
  ((monomer-indices :initarg :monomer-indices :accessor monomer-indices)
   (oligomer-space :initarg :oligomer-space :accessor oligomer-space)))

(cando:make-class-save-load oligomer)

(defun goto-sequence (oligomer index)
  (let* ((bases (loop for monomer across (monomers (oligomer-space oligomer))
                      collect (length (monomers monomer))))
         (digits (sys:positive-integer-to-mixed-base-digits index bases)))
    (loop for monomer-index below (length (monomers (oligomer-space oligomer)))
          for digit in digits
          do (setf (aref (monomer-indices oligomer) monomer-index) digit))))

(defun oligomer-monomer-name-for-monomer (oligomer monomer)
  "Return the monomer name at index in the oligomer"
  (let ((monomer-index (position monomer (monomers oligomer))))
    (elt (monomers monomer) (elt (monomer-indices oligomer) monomer-index))))

(defun oligomer-monomer-name-at-index (oligomer index)
  "Return the monomer name at index in the oligomer"
  (let ((monomer (elt (monomers oligomer) index)))
    (elt (monomers monomer) (elt (monomer-indices oligomer) index))))

(defun oligomer-monomer-names (oligomer)
  "Return a list of monomer names for this oligomer"
  (loop for index below (length (monomers oligomer))
        for monomer across (monomers oligomer)
        collect (elt (monomers monomer) (elt (monomer-indices oligomer) index))))

(defun make-oligomer (oligomer-space &optional (index 0))
  "Build an oligomer in the oligomer space"
  (let ((olig (make-instance 'oligomer
                             :monomer-indices (make-array (length (monomers oligomer-space)) :element-type 'ext:byte32)
                             :oligomer-space oligomer-space)))
    (goto-sequence olig index)
    olig))

(defmethod monomers ((obj oligomer))
  (monomers (oligomer-space obj)))

(defmethod couplings ((obj oligomer))
  (couplings (oligomer-space obj)))


(defun oligomer-space-directional-coupling-iterator-factory (oligomer-space)
  "Return a function that takes an oligomer and returns an
   iterator that iterates over directional couplings in oligomer-space and return (values coupling source-monomer-name target-monomer-name)"
  (let ((all-couplings (loop for coupling across (couplings oligomer-space)
                                   when (typep coupling 'directional-coupling)
                                     collect coupling))
        (monomer-to-index (let ((ht (make-hash-table)))
                            (loop for monomer across (monomers oligomer-space)
                                  for index from 0
                                  do (setf (gethash monomer ht) index))
                            ht)))
    (lambda (oligomer)
      "This is a factory lambda that takes an oligomer and returns an iterator that iterates over the directional couplings in the oligomer"
      (let ((remaining-couplings all-couplings))
        (lambda ()
          "An iterator that iterates over couplings in the oligomer-space/oligomer pair.
              This iterator is closed over an oligomer and it returns the names of monomers in the coupling of the oligomer.
              Couplings are really part of the oligomer-space and so we can precompute them and use them to iterate over multiple oligomers.
              Calling the iterator returns successive (values coupling source-monomer-name target-monomer-name)"
          (when remaining-couplings
            (let* ((coupling (car remaining-couplings))
                   (source-monomer (source-monomer coupling))
                   (source-monomer-index (gethash source-monomer monomer-to-index))
                   (source-monomer-name (elt (monomers source-monomer)
                                             (elt (monomer-indices oligomer) source-monomer-index)))
                   (target-monomer (target-monomer coupling))
                   (target-monomer-index (gethash target-monomer monomer-to-index))
                   (target-monomer-name (elt (monomers target-monomer)
                                             (elt (monomer-indices oligomer) target-monomer-index))))
              (setf remaining-couplings (cdr remaining-couplings))
              (values coupling source-monomer-name target-monomer-name))))))))
  
