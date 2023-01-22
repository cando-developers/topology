(in-package :foldamer)

(define-condition no-matching-context ()
  ((context :initarg :context :accessor context))
  (:report (lambda (condition stream)
             (format stream "no-matching-context ~a" (context condition)))))

(defclass node ()
  ((id :initform (gensym) :initarg :id :accessor id)
   (name :initarg :name :accessor name)
   (spanning-depth :initform nil :initarg :spanning-depth :initform nil :accessor spanning-depth)
   (label :initarg :label :accessor label)))

(cando:make-class-save-load node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t :identity t)
     (format stream "~a :spanning-depth ~a :label ~a" (name obj) (spanning-depth obj) (label obj)))))

(defclass cap-node (node)
  ())

(cando:make-class-save-load cap-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t :identity t)
     (format stream "~a :spanning-depth ~a :label ~a" (name obj) (spanning-depth obj) (label obj)))))

(defclass edge ()
  ((from-node :initarg :from-node :accessor from-node)
   (to-node :initarg :to-node :accessor to-node)
   (name :initarg :name :accessor name)
   (raw-name :initarg :raw-name :accessor raw-name)))

(cando:make-class-save-load edge
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a ~a" (from-node obj) (raw-name obj) (to-node obj)))))

(defun other-node (edge node)
  "If the edge contains node then return the other node, otherwise nil"
  (cond
    ((eq (to-node edge) node)
     (from-node edge))
    ((eq (from-node edge) node)
     (to-node edge))
    (t nil)))

(defclass dag ()
  ((label :initarg :label :accessor label)
   (root :initarg :root :accessor root)
   (nodes :initform nil :accessor nodes)
   (edges :initform nil :accessor edges)))

(cando:make-class-save-load dag
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s" (label obj)))))

(defun recursive-spanning-tree (root depth dag seen)
  (if (gethash root seen)
      nil
      (progn
        (setf (spanning-depth root) depth
              (gethash root seen) t)
        (let ((neighbors (loop for edge in (edges dag)
                               for other-node = (other-node edge root)
                               when other-node
                                 collect other-node)))
          (append neighbors
                  (loop for neighbor in neighbors
                        append (recursive-spanning-tree neighbor (1+ depth) dag seen)))))))

(defun walk-spanning-tree (dag)
  (let ((root (root dag))
        (seen (make-hash-table)))
    (cons root (recursive-spanning-tree root 0 dag seen))))

(defun ensure-monomer (node map)
  (let ((monomer (gethash node map)))
    (unless monomer
      (error "Could not find monomer for node ~a" node))
    monomer))

(defun oligomer-space-from-dag (foldamer dag topology-groups)
  (let ((focus-node (root dag))
        (node-to-monomer (make-hash-table))
        (oligomer-space (make-instance 'topology:oligomer-space
                                       :foldamer foldamer)))
    (loop for node in (nodes dag)
          for name = (name node)
          for names = (gethash name topology-groups)
          for monomer = (make-instance 'topology:monomer
                                       :monomers names
                                       :id name)
          do (setf (gethash node node-to-monomer) monomer)
          do (vector-push-extend monomer (topology:monomers oligomer-space)))
    (loop for edge in (edges dag)
          for source-node = (from-node edge)
          for target-node = (to-node edge)
          for source-monomer = (ensure-monomer source-node node-to-monomer)
          for target-monomer = (ensure-monomer target-node node-to-monomer)
          for coupling-name = (name edge)
          for source-plug-name = (topology:out-plug-name coupling-name)
          for target-plug-name = (topology:in-plug-name coupling-name)
          for directional-coupling = (make-instance 'topology:directional-coupling
                                                    :name coupling-name
                                                    :source-monomer source-monomer
                                                    :target-monomer target-monomer
                                                    :source-plug-name source-plug-name
                                                    :target-plug-name target-plug-name)
          do (setf (gethash source-plug-name (topology:couplings source-monomer))
                   directional-coupling
                   (gethash target-plug-name (topology:couplings target-monomer))
                   directional-coupling)
          do (vector-push-extend directional-coupling (topology:couplings oligomer-space)))
    (values oligomer-space (gethash focus-node node-to-monomer))))

(defun node-from-name (maybe-name label)
  (cond
    ((and (consp maybe-name) (eq (car maybe-name) :cap))
     (make-instance 'cap-node :name (cadr maybe-name) :label label))
    ((symbolp maybe-name)
     (make-instance 'node :name maybe-name :label label))
    (t (error "Illegal name in node-from-name ~s" maybe-name))))

(defun parse-recursive (sub-tree prev-node dag label)
  (cond
    ((null sub-tree))
    ((consp (car sub-tree))
     (parse-recursive (car sub-tree) prev-node dag :head-cons )
     (parse-recursive (cdr sub-tree) prev-node dag :tail-cons))
    ((symbolp (car sub-tree))
     (let* ((plug-name (car sub-tree))
            (node (node-from-name (cadr sub-tree) label)))
       (cond
         ((topology:in-plug-name-p plug-name)
          (let ((edge (make-instance 'edge
                                     :raw-name plug-name
                                     :name (topology:coupling-name plug-name)
                                     :from-node node
                                     :to-node prev-node)))
            (push node (nodes dag))
            (push edge (edges dag))))
         ((topology:out-plug-name-p plug-name)
          (let ((edge (make-instance 'edge
                                     :raw-name plug-name
                                     :name (topology:coupling-name plug-name)
                                     :from-node prev-node
                                     :to-node node)))
            (push node (nodes dag))
            (push edge (edges dag))))
         (t (error "Illegal plug-name ~s" plug-name)))
       (parse-recursive (cddr sub-tree) node dag :symbol)))
    (t (error "Illegal sub-tree ~s" sub-tree))))

(defun parse-label (head)
  head)

(defun parse-dag-for-oligomer-space (labeled-tree)
  (let* ((label (parse-label (car labeled-tree)))
         (tree (cdr labeled-tree))
         (node (node-from-name (car tree) :top))
         (dag (make-instance 'dag :root node :label label)))
    (push node (nodes dag))
    (parse-recursive (cdr tree) node dag :top)
    (walk-spanning-tree dag)
    dag))

(defun validate-dag (dag)
  (let ((label (label dag)))
    (loop for node in (nodes dag)
          for node-name = (name node)
          for depth = (spanning-depth node)
          do (cond
               ((= 1 depth)
                (when (typep node 'cap-node)
                  (warn "dag ~s has a cap-node ~s directly connected to the root" label node-name)))
               ((> depth 1)
                (when (not (typep node 'cap-node))
                  (warn "dag ~s has node ~s at level ~a that should be a cap node or eliminated" label node-name depth)))))))

(defclass training-oligomer-space ()
  ((expression :initarg :expression :accessor expression)
   (expression-dag :initarg :expression-dag :accessor expression-dag)
   (oligomer-space :initarg :oligomer-space :accessor oligomer-space)
   (focus-monomer :initarg :focus-monomer :accessor focus-monomer)
   (monomer-context-matcher :initarg :monomer-context-matcher :accessor monomer-context-matcher)))

(cando:make-class-save-load training-oligomer-space
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s" (expression obj)))))

(defclass foldamer ()
  ((topologys :initform nil :initarg :topologys :accessor topologys)
   (training-oligomer-spaces :initarg :training-oligomer-spaces :accessor training-oligomer-spaces)))

(cando:make-class-save-load foldamer)

(defparameter *foldamers* (make-hash-table))

(defun define-foldamer (name contexts)
  (let ((total-sequences 0)
        (topologys nil)
        (foldamer (make-instance 'foldamer)))
    (let ((training-oligomer-spaces
            (loop for context in contexts
                  collect (let ((dag (parse-dag-for-oligomer-space context)))
                            (validate-dag dag)
                            (multiple-value-bind (oligomer-space focus-monomer)
                                (oligomer-space-from-dag foldamer dag topology::*topology-groups*)
                              #+(or)(format t "trainers ~a ~a~%" (topology:number-of-sequences oligomer-space) context)
                              (let ((topologys-in-oligomer-space (topology:topologys-in-oligomer-space oligomer-space)))
                                (loop for topology in topologys-in-oligomer-space
                                      do (pushnew topology topologys)))
                              (incf total-sequences (topology:number-of-sequences oligomer-space))
                              (let ((monomer-context-matcher (monomer-context:parse context)))
                                (make-instance 'training-oligomer-space
                                               :expression context
                                               :expression-dag dag
                                               :oligomer-space oligomer-space
                                               :focus-monomer focus-monomer
                                               :monomer-context-matcher monomer-context-matcher)))))))
      (reinitialize-instance foldamer
                             :topologys topologys
                             :training-oligomer-spaces training-oligomer-spaces)
      (setf (gethash name *foldamers*)
            foldamer))))


(defun find-oligomer-for-monomer-context (foldamer monomer-context)
  (let ((training-spaces (training-oligomer-spaces foldamer)))
    (loop for training-space in training-spaces
          for oligomer-space = (oligomer-space training-space)
          for focus-monomer = (focus-monomer training-space)
          for monomer-context-matcher = (monomer-context-matcher training-space)
          for num-sequences = (topology:number-of-sequences oligomer-space)
          do (loop for num-seq below num-sequences
                   for oligomer = (topology:make-oligomer oligomer-space num-seq)
                   for match = (monomer-context:match monomer-context-matcher focus-monomer oligomer)
                   for trainer-context = (monomer-context:match-as-string match)
                   when (string= monomer-context trainer-context)
                     do (return-from find-oligomer-for-monomer-context (values oligomer focus-monomer))))))


(defun oligomer-monomer-context-focus-monomer-iterator (training-oligomer-space)
  "Return an iterator that returns successive (values number-remaining oligomer monomer-context focus-monomer focus-monomer-name)
   in the training-oligomer-space. When it runs out it returns nil."
  (let* ((number-of-sequences (topology:number-of-sequences (foldamer:oligomer-space training-oligomer-space)))
         (sequence-index 0)
         (oligomer-space (oligomer-space training-oligomer-space))
         (monomer-context-matcher (monomer-context-matcher training-oligomer-space))
         (focus-monomer (focus-monomer training-oligomer-space))
         (focus-monomer-index (position focus-monomer (topology:monomers oligomer-space))))
    (lambda ()
      (unless (>= sequence-index number-of-sequences)
        (let* ((oligomer (topology:make-oligomer oligomer-space sequence-index))
               (match (monomer-context:match monomer-context-matcher focus-monomer oligomer))
               (monomer-context (monomer-context:match-as-string match))
               (focus-monomer-name (topology:oligomer-monomer-name-at-index oligomer focus-monomer-index)))
          (incf sequence-index)
          (values (- number-of-sequences sequence-index) oligomer monomer-context focus-monomer focus-monomer-name))))))

(defun foldamer-oligomer-monomer-context-focus-monomer-iterator (foldamer)
  "Return an iterator to access all oligomers in the foldamer"
  (let* ((remaining-training-oligomer-spaces (cdr (training-oligomer-spaces foldamer)))
         (cur-training-oligomer-space (car (training-oligomer-spaces foldamer)))
         (inner-iterator (oligomer-monomer-context-focus-monomer-iterator cur-training-oligomer-space)))
    (lambda ()
      (when inner-iterator
        (multiple-value-bind (number-remaining oligomer monomer-context focus-monomer)
            (funcall inner-iterator)
          (cond
            ((> number-remaining 0) (values t oligomer monomer-context focus-monomer))
            ((= number-remaining 0)
             (setf cur-training-oligomer-space (car remaining-training-oligomer-spaces)
                   remaining-training-oligomer-spaces (cdr remaining-training-oligomer-spaces)
                   inner-iterator (and cur-training-oligomer-space (oligomer-monomer-context-focus-monomer-iterator cur-training-oligomer-space)))
             (values t oligomer monomer-context focus-monomer))
            (t nil)))))))

(defun calculate-files (trainer-context &optional root-pathname)
  (let* ((data-dir (if root-pathname (merge-pathnames #P"data/" root-pathname)
                       #P"data/"))
         (output-dir (if root-pathname (merge-pathnames #P"output/" root-pathname)
                         #P"output/"))
         (input-file (make-pathname :name trainer-context :type "input" :defaults data-dir))
         (sdf-file (make-pathname :name trainer-context :type "sdf" :defaults output-dir))
         (internals-file (make-pathname :name trainer-context :type "internals" :defaults output-dir))
         (log-file (make-pathname :name trainer-context :type "log" :defaults output-dir))
         (svg-file (make-pathname :name trainer-context :type "svg" :defaults output-dir))
         (done-file (make-pathname :name trainer-context :type "done" :defaults output-dir)))
    (values input-file done-file sdf-file internals-file log-file svg-file)))

(defun valid-trainer-contexts (foldamer)
  (let ((training-spaces (training-oligomer-spaces foldamer)))
    (loop for training-space in training-spaces
          for oligomer-space = (oligomer-space training-space)
          for focus-monomer = (focus-monomer training-space)
          for monomer-context-matcher = (monomer-context-matcher training-space)
          for num-sequences = (topology:number-of-sequences oligomer-space)
          append (loop for num-seq below num-sequences
                   for oligomer = (topology:make-oligomer oligomer-space num-seq)
                   for match = (monomer-context:match monomer-context-matcher focus-monomer oligomer)
                   unless match
                     do (error "Bad match")
                       collect (monomer-context:match-as-string match)))))

(defun unused-trainer-contexts (foldamer path)
  (let ((valid-trainer-contexts (valid-trainer-contexts foldamer))
        (inputs (directory (merge-pathnames #P"*.input" path))))
    (loop for input in inputs
          for input-context = (pathname-name input)
          unless (find input-context valid-trainer-contexts :test #'string=)
            collect input-context)))

(defun maybe-remove-unused-trainers (foldamer path &key doit)
  (let ((unused-trainer-contexts (unused-trainer-contexts foldamer path)))
    (loop for trainer-context in unused-trainer-contexts
          for files = (multiple-value-list (calculate-files trainer-context))
          do (loop for file in files
                   do (format t "~a ~a~%" (if doit
                                              "Removing"
                                              "Would remove")
                              file)
                   when doit
                   do (delete-file file)))))

(defun generate-training-oligomers (foldamer path &key force-save print)
  (ensure-directories-exist path)
  (let ((foldamer-path (merge-pathnames #P"foldamer.dat" path)))
    (cando:save-cando foldamer foldamer-path))
  (with-open-file (fmakefile (merge-pathnames "makefile" path) :direction :output :if-exists :supersede)
    (let ((all-done-files nil))
      (loop for trainer-context in (valid-trainer-contexts foldamer)
            do (multiple-value-bind (input-file done-file)
                   (calculate-files trainer-context path)
                 (declare (ignore done-file))
                 (multiple-value-bind (local-input-file local-done-file)
                     (calculate-files trainer-context)
                   (ensure-directories-exist input-file)
                   (when (or force-save (null (probe-file input-file)))
                     (with-open-file (fout input-file :direction :output :if-exists :supersede)
                       (format fout "(ql:quickload :topology)~%")
                       (format fout "(format t \"Building trainer in file ~~s~~%\" *load-pathname*)~%")
                       (format fout "(defparameter agg (foldamer:build-trainer ~s))~%" trainer-context)
                       (format fout "(sys:exit 0)~%")
                       (when print
                         (format t "Generating trainer for ~a~%" trainer-context))
                       )
                     (push local-done-file all-done-files)
                     (format fmakefile "~a : ~a~%" local-done-file local-input-file)
                     (format fmakefile "~a$(CLASP) -t c -f cclasp -l ~s~%~%" #\tab (namestring local-input-file))
                     ))))
      (format fmakefile "all: ~{~a ~}~%" all-done-files)
      )))

(defun register-topologys (foldamer)
  (loop for topology in (topologys foldamer)
        do (cando:register-topology topology (topology:name topology))))

(defun extract-one-internal (joint flog)
  (let ((name (kin:joint/name joint)))
    (cond
      ((typep joint 'kin:jump-joint)
       (make-instance 'topology:jump-internal
                      :name name))
      ((typep joint 'kin:complex-bonded-joint)
       (let ((distance (kin:bonded-joint/get-distance joint))
             (angle (kin:bonded-joint/get-theta joint))
             (dihedral (kin:bonded-joint/get-phi joint)))
         (make-instance 'topology:complex-bonded-internal
                        :name name
                        :bond distance
                        :angle angle
                        :dihedral dihedral)))
      ((typep joint 'kin:bonded-joint)
       (let ((distance (kin:bonded-joint/get-distance joint))
             (angle (kin:bonded-joint/get-theta joint))
             (dihedral (kin:bonded-joint/get-phi joint)))
         (make-instance 'topology:bonded-internal
                        :name name
                        :bond distance
                        :angle angle
                        :dihedral dihedral)))
      (t (format flog "unknown-joint ~a ~a~%" name joint)))))


(defun out-of-focus-atresidue-internals (atmolecule focus-atresidue flog)
  "Extract the out-of-focus internals, these are internal coordinates in atresidues that are connected to the focus-atresidue
by through a parent or grandparent"
  (let ((joint-to-atresidue (make-hash-table)))
    (loop with atresidues = (topology:atresidues atmolecule)
          for atresidue-index below (length atresidues)
          for atresidue = (aref atresidues atresidue-index)
          do (loop with joints = (topology:joints atresidue)
                   for joint-index below (length joints)
                   for joint = (aref joints joint-index)
                   do (setf (gethash joint joint-to-atresidue) atresidue)))
    (let ((out-of-focus-internals (make-hash-table :test 'equal)))
      (loop with atresidues = (topology:atresidues atmolecule)
            for atresidue-index below (length atresidues)
            for atresidue = (aref atresidues atresidue-index)
            do (unless (eq atresidue focus-atresidue) ; skip focus-atresidue
                 (let* ((joints (topology:joints atresidue))
                        (first-joint (elt joints 0))
                        (first-parent (kin:get-parent first-joint))
                        (first-parent-atresidue (gethash first-parent joint-to-atresidue)))
                   (when (eq first-parent-atresidue focus-atresidue) ; if true then atresidue follows focus-atresidue
                     (let ((internals-vector (loop
                                               named out-of-focus-internals
                                               with internals-vector = (make-array 8 :adjustable t :fill-pointer 0)
                                               for joint-index below (length joints)
                                               for joint = (aref joints joint-index)
                                               do (let* ((parent (kin:get-parent joint))
                                                         (grandparent (when parent (kin:get-parent parent)))
                                                         (great-grandparent (when grandparent (kin:get-parent grandparent)))
                                                         (parent-atresidue (and parent (gethash parent joint-to-atresidue)))
                                                         (grandparent-atresidue (and grandparent (gethash grandparent joint-to-atresidue)))
                                                         (great-grandparent-atresidue (and great-grandparent (gethash great-grandparent joint-to-atresidue))))
                                                    (if (or (eq focus-atresidue parent-atresidue)
                                                            (eq focus-atresidue grandparent-atresidue)
                                                            #+(or)(eq focus-atresidue great-grandparent-atresidue))
                                                        (vector-push-extend (extract-one-internal joint flog) internals-vector)
                                                        (return-from out-of-focus-internals internals-vector))))))
                       (setf (gethash (topology:stereoisomer-name atresidue) out-of-focus-internals) internals-vector))))))
      out-of-focus-internals)))

(defun extract-focus-atresidue-internals (focus-atresidue atmolecule total-count flog)
  "Extract the internal coordinates for the atresidue.
Also extract the out-of-focus-internals, the internals that leave the focus residue but have a parent or grandparent in the focus residue.
We need these to match fragment internals with each other later."
  (let* ((internals (loop for joint across (topology:joints focus-atresidue)
                          collect (extract-one-internal joint flog)))
         (out-of-focus-internals (out-of-focus-atresidue-internals atmolecule focus-atresidue flog)))
    (make-instance 'topology:fragment-internals
                   :index total-count
                   :internals (coerce internals 'vector)
                   :out-of-focus-internals out-of-focus-internals)))


(defun build-trainer (foldamer trainer-context &key (steps 3) (load-pathname *load-pathname*) build-info-msg verbose)
  (let ((root-pathname (make-pathname :directory (butlast (pathname-directory load-pathname)))))
    (multiple-value-bind (input-file done-file sdf-file internals-file log-file svg-file)
        (calculate-files trainer-context root-pathname)
      (declare (ignore input-file))
      (ensure-directories-exist sdf-file)
      (with-open-file (fsdf sdf-file :direction :output
                                     :if-exists :append
                                     :if-does-not-exist :create)
        (with-open-file (flog log-file :direction :output
                                       :if-exists :append
                                       :if-does-not-exist :create)
          (ensure-directories-exist log-file)
          (multiple-value-bind (oligomer focus-monomer)
              (find-oligomer-for-monomer-context foldamer trainer-context)
            (let* ((fragment-conformations (if (probe-file internals-file)
                                               (topology:load-fragment-conformations internals-file)
                                               (make-instance 'topology:fragment-conformations
                                                              :focus-monomer-name (topology:current-stereoisomer-name focus-monomer oligomer)
                                                              :monomer-context trainer-context))))
              (register-topologys foldamer)
              #+(or)(format flog "Building trainer for monomer context ~a~%" (dump-local-monomer-context focus-monomer))
              (let* ((conf (topology:make-conformation oligomer :focus-monomer focus-monomer))
                     (agg (topology:aggregate conf))
                     (molecule (cando:mol agg 0))
                     (number-of-atoms (chem:number-of-atoms molecule))
                     (monomer-positions (topology:monomer-positions conf))
                     (focus-residue-index (gethash focus-monomer monomer-positions))
                     (ataggregate (topology:ataggregate conf))
                     (atmolecule (aref (topology:atmolecules ataggregate) 0))
                     (focus-atresidue (aref (topology:atresidues atmolecule) focus-residue-index))
                     (total-count (topology:total-count fragment-conformations)))
                (with-open-file (fout svg-file :direction :output :if-exists :supersede)
                  (let* ((sketch2d (sketch2d:sketch2d molecule))
                         (svg (sketch2d:svg sketch2d)))
                    (cl-svg:stream-out fout (sketch2d:render-svg-scene svg))))
                (when (plusp (- steps total-count))
                  (format flog "internals ~a~%" trainer-context)
                  (loop for count below (- steps total-count)
                        do (progn
                             (block once
                               (handler-bind
                                   ((chem:minimizer-error (lambda (err)
                                                            (let ((save-filename (make-pathname :name (format nil "~a-~a" (pathname-name flog) count)
                                                                                                :type "cando"
                                                                                                :defaults flog)))
                                                              (format flog "build-trainer - the minimizer reported: ~a - writing to ~a~%" err save-filename)
                                                              (when verbose (format t "Encountered a minimizer error: ~a~%" err))
                                                              (invoke-restart 'cando:save-and-skip-rest-of-minimization save-filename)
                                                              (return-from once nil))))
                                    (smirnoff:missing-dihedral (lambda (err)
                                                                 (let ((save-filename (make-pathname :type "cando" :defaults flog)))
                                                                   (format flog "Missing dihedral ~a - saving molecule to ~s~%" err save-filename)
                                                                   (cando:save-cando (smirnoff:molecule err) save-filename))
                                                                 (signal err))))
                                 (progn
                                   (cando:starting-geometry-with-restarts agg :verbose verbose)))
                               (format flog "Found a starting geometry for total-count: ~a~%" total-count)
                               (let ((data-items nil)
                                     (maybe-bad-geometry (topology:bad-geometry-p agg)))
                                 (if (null maybe-bad-geometry)
                                     (progn
                                       (format flog "Passed (not (topology:bad-geometry-p agg))~%")
                                       (push (cons "status" "good-geometry") data-items)
                                       (topology::copy-atom-positions-into-joints conf)
                                       (topology::update-joint-tree-internal-coordinates conf)
                                       (let* ((fragment-internals (extract-focus-atresidue-internals focus-atresidue atmolecule total-count flog))
                                              (seen-index (topology:seen-fragment-internals fragment-conformations fragment-internals)))
                                         (if (topology:good-fragment-internals fragment-internals)
                                             (if (not seen-index)
                                                 (progn
                                                   (push (cons "Seen" "false") data-items)
                                                   (push (cons "internals-index" (length (topology:fragments fragment-conformations))) data-items)
                                                   (format flog "Saving fragment internals for conformation: ~a~%" total-count)
                                                   (topology:dump-fragment-internals fragment-internals flog)
                                                   (push fragment-internals (topology:fragments fragment-conformations)))
                                                 (progn
                                                   (push (cons "Seen" "true") data-items)
                                                   (format flog "Ignoring conformation ~a - seen before at ~a~%" total-count seen-index)
                                                   (topology:dump-fragment-internals fragment-internals flog)
                                                   ))
                                             (progn
                                               (push (cons "status" "ignoring conformation - failed (topology:good-fragment-internals fragment-internals)") data-items)
                                               (format flog "Ignoring conformation ~a - failed (topology:good-fragment-internals fragment-internals)~%" total-count)
                                               (topology:dump-fragment-internals fragment-internals flog)
                                               ))))
                                     (progn
                                       (push (cons "status" (format nil "Failed (not (topology:bad-geometry-p agg)) for conformation: ~a~%problem: ~a~%" total-count maybe-bad-geometry)) data-items)
                                       (format flog "Failed (not (topology:bad-geometry-p agg)) for conformation: ~a~%problem: ~a~%" total-count maybe-bad-geometry)))
                                 (sdf:write-sdf-stream agg fsdf :name trainer-context :data-items (list* (cons "total-count index" total-count) data-items)))
                               (incf total-count)))))
                (setf (topology:total-count fragment-conformations) total-count)
                (topology:save-fragment-conformations fragment-conformations internals-file)
                (with-open-file (fout done-file :direction :output :if-exists :supersede)
                  (format fout "(:number-of-atoms ~a " number-of-atoms)
                  (when build-info-msg
                    (format fout ":build-info-msg ~s " build-info-msg))
                  (format fout ":finished-steps ~a " total-count)
                  (let ((fragment-conformations-length (length (topology:fragments fragment-conformations))))
                    (format fout ":hits ~a " fragment-conformations-length)
                    (if (/= total-count 0)
                        (format fout ":hit-to-total-ratio ~5,4f " (/ (float fragment-conformations-length 1.0d0) (float total-count 1.0d0)))
                        (format fout ":hit-to-total-ratio 0.0 "))
                    (format fout ")~%")))))))))))


(defun prepare-to-build-trainer (&key (smirnoff #P"~/Development/openff-sage/inputs-and-results/optimizations/vdw-v1/forcefield/force-field.offxml" ))
  (leap:load-smirnoff-params smirnoff)
  )

(defun extract-fragment-conformations-map (filename &key (parallel nil))
  "Parallel version is slower than serial because of false sharing in GC"
  (let* ((foldamer-filename (merge-pathnames #P"foldamer.dat" filename))
         (foldamer (cando:load-cando foldamer-filename))
         (fragment-conformations-map (make-instance 'topology:fragment-conformations-map)))
    (flet ((extract-one (trainer-name)
             (multiple-value-bind (input-file done-file sdf-file internals-file)
                 (calculate-files trainer-name filename)
               (declare (ignore input-file done-file sdf-file))
               (if (probe-file internals-file)
                   (let ((fragment-conformations (topology:load-fragment-conformations internals-file)))
                     (cons trainer-name fragment-conformations))
                   (warn "Could not read file ~a" internals-file)))))
      (let ((names-values
              (if parallel
                  (lparallel:pmapcar #'extract-one (valid-trainer-contexts foldamer))
                  (mapcar #'extract-one (valid-trainer-contexts foldamer)))))
        (loop for name-result in names-values
              for trainer-name = (car name-result)
              for fragment-conformations = (cdr name-result)
              do (when trainer-name
                   (setf (gethash trainer-name (topology:monomer-context-to-fragment-conformations fragment-conformations-map))
                         fragment-conformations)))
        fragment-conformations-map))))


(defun check-fragment-conformations-map (filename)
  (let* ((foldamer-filename (merge-pathnames #P"foldamer.dat" filename))
         (foldamer (cando:load-cando foldamer-filename))
         (fragment-conformations-map (make-instance 'topology:fragment-conformations-map))
	 (left-count 0)
	 (right-count 0))
    (loop for trainer-name in (valid-trainer-contexts foldamer)
	  do (multiple-value-bind (input-file done-file sdf-file internals-file)
				  (calculate-files trainer-name filename)
				  (declare (ignore input-file done-file sdf-file))
				  (if (probe-file internals-file)
				      (let ((fragment-conformations (topology:load-fragment-conformations internals-file)))
					(loop for fragment-internals in (topology:fragments fragment-conformations)
					      for internals = (topology:internals fragment-internals)
					      for internal0 = (first internals)
					      when (eq (topology:name internal0) :cm)
					      do (let* ((j1 (second internals))
							(j2 (third internals))
							(j3 (fourth internals))
							(d1 (/ (topology:dihedral j1) 0.0174533))
							(d2 (/ (topology:dihedral j2) 0.0174533))
							(d3 (/ (topology:dihedral j3) 0.0174533))
							(ad12 (foldamer::angle-difference d1 d2))
							(ad23 (foldamer::angle-difference d2 d3))
							(ad31 (foldamer::angle-difference d3 d1)))
						   (if (> ad12 0)
						       (incf right-count)
						     (incf left-count))
						   (format t "~a ->  cm ~a ~a ~a : ~6,1f ~6,1f ~6,1f~%"
							   trainer-name
							   (topology:name j1)
							   (topology:name j2)
							   (topology:name j3)
							   ad12
							   ad23
							   ad31)))))))
    (format t " Left turning count: ~a~%" left-count)
    (format t "Right turning count: ~a~%" right-count)))


(defun recursive-dump-local-monomer-context (monomer prev-monomer depth)
  (when (>= depth 0)
    (let (rest)
      (maphash (lambda (key coupling)
                 (let ((other-monomer (topology:other-monomer coupling monomer)))
                   (unless (eq other-monomer prev-monomer)
                     (let ((other-tree (recursive-dump-local-monomer-context other-monomer monomer (1- depth))))
                       (when other-tree
                         (push (cons key other-tree) rest))))
                   ))
                 (topology:couplings monomer))
               (let ((sorted-rest (sort rest #'string< :key (lambda (val) (string (car val))))))
                 (list* (topology:monomers monomer) (mapcar (lambda (alist) (list (car alist) (cdr alist))) sorted-rest))))))

(defun dump-local-monomer-context (monomer)
  (recursive-dump-local-monomer-context monomer nil 2))

(defmethod topology:foldamer-monomer-context (monomer oligomer-or-space (foldamer foldamer))
  "Return the monomer-context string for the monomer in the oligomer that is made of the foldamer"
  (unless (find monomer (topology:monomers oligomer-or-space))
    (error "The monomer ~a was not found in the oligomer-or-space ~a" monomer oligomer-or-space))
  (block monomer-context
    (loop for training-oligomer-space in (training-oligomer-spaces foldamer)
          for monomer-context-matcher = (monomer-context-matcher training-oligomer-space)
          for match = (monomer-context:match monomer-context-matcher monomer oligomer-or-space)
          when match
            do (progn
                 (return-from monomer-context (monomer-context:match-as-string match))))
    (error 'no-matching-context
           :context (dump-local-monomer-context monomer))))


(defun verify-foldamer-describes-oligomer-space (foldamer oligomer-space &key print)
  "Check that every monomer in the oligomer space has a monomer-context within the foldamer"
  (let ((used-contexts-set (make-hash-table :test 'equal)))
    (loop for monomer across (topology:monomers oligomer-space)
          for monomer-context = (topology:foldamer-monomer-context monomer oligomer-space foldamer)
          do (format t "monomer-context = ~s~%" monomer-context)
          do (setf (gethash monomer-context used-contexts-set) t)
          if monomer-context
            do (when print
                 (format t "monomer-context ~a~%   is matched by ~a~%" (recursive-dump-local-monomer-context monomer nil 1) monomer-context))
          else
            do (error "Foldamer does not describe oligomer space"))
    (let (used-contexts unused-contexts)
      (maphash (lambda (key value)
                 (declare (ignore value))
                 (if (gethash key used-contexts-set)
                     (push key used-contexts)
                     (push key unused-contexts)))
               used-contexts-set)
      (values t used-contexts unused-contexts))))
