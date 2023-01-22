(in-package :foldamer)



(defun done-data (trainer-context)
  (multiple-value-bind (input-file done-file)
    (foldamer:calculate-files trainer-context)
    (declare (ignore input-file))
    (if (null (probe-file done-file))
	(values nil)
      (with-open-file (done-fin done-file :direction :input)
        (read done-fin)))))

(defun needs-work (trainer-context steps)
  "Return (values needs-work hits) 
    If this trainer-context has not run enough steps then needs-work is T.
    Also return the number of hits that have been found."
  (let* ((done-data (done-data trainer-context))
	 (done-steps (getf done-data :finished-steps))
	 (hits (getf done-data :hits)))
    (or (null done-data) (or (< done-steps steps) (= hits 0)))))

(defun max-finished-steps (spiros)
  (let ((max-finished-steps 0)
        (valid-trainer-contexts (foldamer:valid-trainer-contexts spiros)))
    (loop for context in valid-trainer-contexts
          do (multiple-value-bind (input-file done-file)
                 (foldamer:calculate-files context)
               (declare (ignore input-file))
               #+(or)(format t "Reading ~a~%" done-file)
               (when (probe-file done-file)
                 (with-open-file (fin done-file)
                   (let ((info (read fin)))
                     (setf max-finished-steps (max max-finished-steps (getf info :finished-steps)))
                     )))))
    max-finished-steps))


(defun needs-extract-conformations (path spiros)
  (let ((valid-trainer-contexts (foldamer:valid-trainer-contexts spiros))
        (conformations-pathname  path))
    (if (null (probe-file conformations-pathname))
        t ;; Need to update conformations if file doesn't exist
        (let ((latest-done-file-write-date 0)
              (conformations-write-date (file-write-date conformations-pathname)))
          (loop for context in valid-trainer-contexts
                do (multiple-value-bind (input-file done-file)
                       (foldamer:calculate-files context)
                     (declare (ignore input-file))
                     #+(or)(format t "Reading ~a~%" done-file)
                     (when (probe-file done-file)
                       (setf latest-done-file-write-date (max latest-done-file-write-date (file-write-date done-file))))))
          (< conformations-write-date latest-done-file-write-date)))))

(defun foldamer-update-trainers (spiros &key remove-unused force-save)
  (foldamer:maybe-remove-unused-trainers spiros #P"./" :doit remove-unused)
  (foldamer:generate-training-oligomers spiros #P"./" :force-save force-save))

(defun foldamer-status (&key spiros)
  (let* ((trainer-contexts (foldamer:valid-trainer-contexts spiros))
	 (trainer-done-data (make-hash-table :test 'equal))
	 (max-finished-steps (max-finished-steps spiros)))
    (loop for trainer-context in trainer-contexts
	  for done-data = (done-data trainer-context)
	  do (setf (gethash trainer-context trainer-done-data) done-data))
    (loop for trainer-context in trainer-contexts
	  for done-data = (gethash trainer-context trainer-done-data)
	  for hits = (getf done-data :hits)
	  for finished-steps = (getf done-data :finished-steps)
	  when (and hits (= hits 0))
            do (format t "~a :finished-steps ~a :hits ~a~%" trainer-context finished-steps hits))
    (let ((number-of-trainers-with-max-finished-steps 
	   (loop for trainer-context in trainer-contexts
		 for done-data = (gethash trainer-context trainer-done-data)
		 for finished-steps = (getf done-data :finished-steps)
		 when (and done-data (= finished-steps max-finished-steps))
		   count 1)))
      (format t "Number of trainers with ~a finished steps = ~a~%" 
	      max-finished-steps
	      number-of-trainers-with-max-finished-steps))
    (format t "Number of trainers ~a~%" (length trainer-contexts))
    ))

(defclass trainer-job ()
  ((trainer-context :initarg :trainer-context :reader trainer-context)
   (number-of-atoms :initarg :number-of-atoms :reader number-of-atoms)
   (node-index :initarg :node-index :reader node-index)
   ))

(cando:make-class-save-load
 trainer-job
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (trainer-context obj) (number-of-atoms obj)))))

(defun foldamer-setup (spiros
		       &key
		       (nodes 12)
		       (clasp "$CLASP -l code.lisp")
		       (script "jobs.dat")
		       (finished-steps 2 finished-steps-p)
		       (advance-steps 10 advance-steps-p))
  (foldamer-update-trainers spiros :remove-unused t)
  (let* ((max-finished-steps (max-finished-steps spiros))
         (steps (cond
		 ((and (null advance-steps-p) finished-steps-p) finished-steps)
		 ((and (null finished-steps-p) advance-steps-p) (+ max-finished-steps advance-steps))
		 ((and (null finished-steps-p) (null advance-steps-p))
		  max-finished-steps)
		 (t (error "You must provide only one option advance-steps or finished-steps or neither"))))
         (trainer-contexts (foldamer:valid-trainer-contexts spiros))
         (node-total nodes)
	 (jobs-per-node (ceiling (length trainer-contexts) node-total))
         (node-work (make-hash-table :test 'eql))
         (trainer-count 0)
         (unsorted-trainers-that-need-work (loop for trainer-context in trainer-contexts
						 when (needs-work trainer-context steps)
						 collect (let* ((oligomer (find-oligomer-for-monomer-context spiros trainer-context))
								(molecule (topology:build-molecule oligomer))
								(number-of-atoms (chem:number-of-atoms molecule)))
							   (cons number-of-atoms trainer-context))))
         (trainers-that-need-work (sort unsorted-trainers-that-need-work #'< :key #'car)))
    ;; Note: Sort order is reverse of what we want because below we push the trainer-context into list
    (format t "Maximum finished-steps ~a~%" max-finished-steps)
    (format t "New steps              ~a~%" steps)
    (loop for pair in trainers-that-need-work
          for number-of-atoms = (car pair)
          for trainer-file = (cdr pair)
          for index from 0
          for node-index = (mod index node-total)
          for trainer-context = (pathname-name (pathname trainer-file))
          do (progn
               #+(or)(format t "number-of-atoms ~a node-index ~a~%" number-of-atoms node-index)
               (let ((trainer-job (make-instance 'trainer-job
                                                 :trainer-context trainer-context
                                                 :number-of-atoms number-of-atoms
                                                 :node-index node-index)))
                 (push trainer-job (gethash node-index node-work nil))
                 (format t "trainer-job ~a~%" trainer-job))
               (incf trainer-count)
               (when (gethash node-index node-work)
                 #+(or)(format t "~a ~a~%" node-index trainer-context))))
    (with-open-file (fout script :direction :output :if-exists :supersede)
		    (loop for node-index below (min jobs-per-node (hash-table-count node-work))
			  do (format fout "~a -e \"(foldamer:foldamer-run-node ~a :steps ~a)\"~%"
				     clasp
				     node-index
				     steps)))
    (cando:save-cando node-work "node-work.cando")
    (format t "~a trainers need to be trained~%" (if (= trainer-count 0) "No" trainer-count))
    )
  (sys:quit))


(defun foldamer-extract-conformations (&key (path #P"./conformations.cando") spiros verbose)
  (format t "Extracting conformations for the path: ~a~%" path)
  (let ((foldamer-conformations-map (foldamer:extract-fragment-conformations-map path)))
    (format t "Matching fragment conformations~%")
    (let ((matched (foldamer:optimize-fragment-conformations-map foldamer-conformations-map spiros verbose)))
      (multiple-value-bind (fragment-total matched-fragments-total missing-matches-total)
          (topology:matched-fragment-conformations-summary matched)
        (format t "Total fragments:    ~8d~%" fragment-total)
        (format t "Matching fragments: ~8d~%" matched-fragments-total)
        (format t "Missing matches:    ~8d~%" missing-matches-total)
        (format t "Fraction missing:   ~8,5f~%" (float (/ missing-matches-total (+ missing-matches-total matched-fragments-total)))))
      (format t "Saving ~a~%" path)
      (cando:save-cando matched path)
      )))


(defun foldamer-describe-missing-fragment-matches (fragment-conformations-map)
  (let ((*print-pretty* nil))
    (maphash (lambda (key value)
               (let ((before-monomer-context (aref (topology:monomer-contexts-vector fragment-conformations-map) (topology:fragment-match-key-before-monomer-context-index key)))
                     (after-monomer-context (aref (topology:monomer-contexts-vector fragment-conformations-map) (topology:fragment-match-key-after-monomer-context-index key))))
                 (format t "~a ~a ~a~%" before-monomer-context after-monomer-context value)))
             (topology:missing-fragment-matches fragment-conformations-map))))

(defun dump-internals (message internals &optional compare)
  (format t "~a~%" message)
  (loop for internal in (coerce internals 'list)
        for index below (if compare (length compare) (length internals))
        for name = (topology:name internal)
        for dih = (/ (topology:dihedral internal) 0.0174533)
        do (progn
             (format t " ~4a ~7,1f" name dih)
             (when compare
               (let ((other-dih (/ (topology:dihedral (elt compare index)) 0.0174533)))
                 (format t "   ~7,1f  ~a" other-dih (if (< (foldamer:angle-difference dih other-dih) 30.0) "MATCH" "---"))))
             (format t "~%"))))
             

(defun foldamer-describe-missing-match (fragment-conformations-map before-monomer-context after-monomer-context &optional before-fragment-index after-monomer)
  (let* ((before-monomer-context-index (gethash before-monomer-context (topology:monomer-context-index-map fragment-conformations-map)))
         (after-monomer-context-index (gethash after-monomer-context (topology:monomer-context-index-map fragment-conformations-map)))
         (match-key (topology:make-fragment-match-key
                     :before-monomer-context-index before-monomer-context-index
                     :after-monomer-context-index after-monomer-context-index))
         (missing-match-vec (gethash match-key (topology:missing-fragment-matches fragment-conformations-map))))
    (cond
      ((null before-fragment-index)
       (format t "vec ~a~%" missing-match-vec))
      ((null after-monomer)
       (let* ((fragment-conformations (gethash before-monomer-context (topology:monomer-context-to-fragment-conformations fragment-conformations-map)))
              (fragment-conformation (elt (topology:fragments fragment-conformations) before-fragment-index))
              (out-of-focus (topology:out-of-focus-internals fragment-conformation)))
         (maphash (lambda (key value)
                    (format t "following monomer ~s ~s~%" key value))
                  out-of-focus)))
      (t 
       (let* ((fragment-conformations (gethash before-monomer-context (topology:monomer-context-to-fragment-conformations fragment-conformations-map)))
              (fragment-conformation (elt (topology:fragments fragment-conformations) before-fragment-index))
              (out-of-focus (topology:out-of-focus-internals fragment-conformation))
              (out-of-focus-internals (gethash after-monomer out-of-focus)))
         (if out-of-focus-internals
             (progn
               (dump-internals "before internals" out-of-focus-internals)
               (loop with after-fragment-conformations = (gethash after-monomer-context (topology:monomer-context-to-fragment-conformations fragment-conformations-map))
                     for after-fragment-index below (length (topology:fragments after-fragment-conformations))
                     for after-fragment-conformation = (elt (topology:fragments after-fragment-conformations) after-fragment-index)
                     for after-internals = (topology:internals after-fragment-conformation)
                     do (dump-internals (format nil "after internals #~a" after-fragment-index) after-internals out-of-focus-internals)))))))))

(defun foldamer-check-conformations (&optional (path #P"./") spiros)
  (format t "checking conformations~%")
  (let ((foldamer-conformations-map (foldamer::check-fragment-conformations-map path)))
    ))

(defun parallel-map (fn jobs)
  (let* ((channel (lparallel:make-channel))
         (tasks (loop for job in jobs
                      collect (lparallel:submit-task channel
                                                     (let ((*job* job))
                                                      (lambda ()
                                                        (funcall fn *job*)))))))
    (loop for task in tasks
          do (lparallel:receive-result channel))))

(defun foldamer-run-node (node-index &key (steps 2) testing (parallel t) (no-fail t))
  (format t "Loading force-field~%")
  (foldamer:prepare-to-build-trainer :smirnoff "./force-field.offxml")
  (let* ((all-node-work (cando:load-cando "node-work.cando"))
         (node-trainers (gethash node-index all-node-work))
         (foldamer-dat-pathname (merge-pathnames #P"./foldamer.dat"))
         (foldamer (cando:load-cando foldamer-dat-pathname)))
    (flet ((one-trainer (trainer-job node-index)
             (let* ((trainer-context (trainer-context trainer-job))
                    (input-file (merge-pathnames (make-pathname :directory '(:relative "data")
                                                                :name trainer-context
                                                                :type "input")))
                    (worker-index (lparallel:kernel-worker-index))
                    (worker-log-name (make-pathname :name (format nil "worker-~2,'0d-~2,'0d" node-index worker-index)
                                                    :directory (list :relative "worker-logs"))))
               (ensure-directories-exist worker-log-name)
               (with-open-file (worker-log worker-log-name :direction :output :if-exists :append :if-does-not-exist :create)
                 (let ((*standard-output* worker-log))
                   (format t "About to build-trainer: ~a~%" trainer-job)
                   (if no-fail
                       (handler-bind
                           ((error (lambda (err)
                                     (format t "~s~%" err)
                                     (let ((clasp-debug:*frame-filters* nil))
                                       (clasp-debug:print-backtrace)))))
                         (foldamer:build-trainer foldamer trainer-context :load-pathname input-file :steps steps :build-info-msg (list :node-index node-index :verbose nil)))
                       (foldamer:build-trainer foldamer trainer-context :load-pathname input-file :steps steps :verbose nil)))))))
      (if parallel
          (parallel-map (lambda (trainer) (one-trainer trainer node-index)) node-trainers)
          (mapc (lambda (trainer) (one-trainer trainer node-index)) node-trainers))))
  (unless testing (sys:quit)))

