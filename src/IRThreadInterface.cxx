// -*- Mode: c++; tab-width: 8; c-basic-offset: 2; indent-tabs-mode: t -*-
// NOTE: the first line of this file sets up source code indentation rules
// for Emacs; it is also a hint to anyone modifying this file.

/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


// File         : the_thread_interface.cxx
// Author       : Pavel A. Koshevoy
// Created      : Fri Feb 16 10:11:00 MST 2007
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : An abstract thread class interface.

// local includes:
#include "IRThreadInterface.h"
#include "IRThreadPool.h"
#include "IRMutexInterface.h"
#include "IRTerminator.h"
#include "IRTransaction.h"
#include "itkIRUtils.h"

// system includes:
#include <iostream>

// namespace access:
using std::cout;
using std::cerr;
using std::endl;

//----------------------------------------------------------------
// DEBUG_THREAD
// 
// #define DEBUG_THREAD


//----------------------------------------------------------------
// MUTEX
// 
static the_mutex_interface_t * MUTEX()
{
  static the_mutex_interface_t * mutex_ = nullptr;
  if (mutex_ == nullptr)
  {
    mutex_ = the_mutex_interface_t::create();
  }
  
  return mutex_;
}


//----------------------------------------------------------------
// the_thread_interface_t::creator_
// 
the_thread_interface_t::creator_t
the_thread_interface_t::creator_ = nullptr;

//----------------------------------------------------------------
// the_thread_interface_t::the_thread_interface_t
// 
the_thread_interface_t::the_thread_interface_t(the_mutex_interface_t * mutex):
  mutex_(mutex),
  stopped_(true),
  sleep_when_idle_(false),
  sleep_microsec_(10000),
  active_transaction_(nullptr),
  thread_pool_(nullptr),
  thread_pool_cb_data_(nullptr)
{
  the_lock_t<the_mutex_interface_t> locker(MUTEX());
  static unsigned int counter = 0;
  id_ = counter++;
}

//----------------------------------------------------------------
// the_thread_interface_t::~the_thread_interface_t
// 
the_thread_interface_t::~the_thread_interface_t()
{
  mutex_->delete_this();
}

//----------------------------------------------------------------
// the_thread_interface_t::set_creator
// 
void
the_thread_interface_t::set_creator(the_thread_interface_t::creator_t creator)
{
  creator_ = creator;
}

//----------------------------------------------------------------
// the_thread_interface_t::create
// 
the_thread_interface_t *
the_thread_interface_t::create()
{
  if (creator_ == nullptr) return nullptr;
  return creator_();
}

//----------------------------------------------------------------
// the_thread_interface_t::set_mutex
// 
void
the_thread_interface_t::set_mutex(the_mutex_interface_t * mutex)
{
  mutex_->delete_this();
  mutex_ = mutex;
}

//----------------------------------------------------------------
// the_thread_interface_t::set_idle_sleep_duration
// 
void
the_thread_interface_t::set_idle_sleep_duration(bool enable,
						unsigned int microseconds)
{
  sleep_when_idle_ = enable;
  sleep_microsec_ = microseconds;
}

//----------------------------------------------------------------
// the_thread_interface_t::push_back
// 
void
the_thread_interface_t::push_back(the_transaction_t * transaction)
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  transactions_.push_back(transaction);
}

//----------------------------------------------------------------
// the_thread_interface_t::push_front
// 
void
the_thread_interface_t::push_front(the_transaction_t * transaction)
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  transactions_.push_front(transaction);
}

//----------------------------------------------------------------
// the_thread_interface_t::push_back
// 
void
the_thread_interface_t::push_back(std::list<the_transaction_t *> & schedule)
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  transactions_.splice(transactions_.end(), schedule);
}

//----------------------------------------------------------------
// the_thread_interface_t::has_work
// 
bool
the_thread_interface_t::has_work() const
{
  return (active_transaction_ != nullptr) || !transactions_.empty();
}

//----------------------------------------------------------------
// the_thread_interface_t::start
// 
void
the_thread_interface_t::start(the_transaction_t * transaction)
{
  push_back(transaction);
  start();
}

//----------------------------------------------------------------
// the_thread_interface_t::stop
// 
void
the_thread_interface_t::stop()
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  if (!stopped_)
  {
#ifdef DEBUG_THREAD
    cerr << "stopping thread: " << this << endl;
#endif
    stopped_ = true;
    terminators().terminate();
  }
}

//----------------------------------------------------------------
// the_thread_interface_t::flush
// 
void
the_thread_interface_t::flush()
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  while (!transactions_.empty())
  {
    the_transaction_t * t = remove_head(transactions_);
#ifdef DEBUG_THREAD
    cerr << this << " flush " << t << endl;
#endif
    t->notify(this, the_transaction_t::SKIPPED_E);
  }
}

//----------------------------------------------------------------
// the_thread_interface_t::terminate_transactions
// 
void
the_thread_interface_t::terminate_transactions()
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  
  // remove any further pending transactions:
  while (!transactions_.empty())
  {
    the_transaction_t * t = remove_head(transactions_);
#ifdef DEBUG_THREAD
    cerr << this << " terminate_transactions " << t << endl;
#endif
    t->notify(this, the_transaction_t::SKIPPED_E);
  }
  
  terminators().terminate();
}

//----------------------------------------------------------------
// the_thread_interface_t::stop_and_go
// 
void
the_thread_interface_t::stop_and_go(the_transaction_t * transaction)
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  
  // remove any further pending transactions:
  while (!transactions_.empty())
  {
    the_transaction_t * t = remove_head(transactions_);
#ifdef DEBUG_THREAD
    cerr << this << " stop_and_go (1) " << t << endl;
#endif
    t->notify(this, the_transaction_t::SKIPPED_E);
  }
  
  // terminate the currently executing transaction:
  terminators().terminate();
  
  // schedule the next transaction:
  transactions_.push_back(transaction);
  
  // start the thread if it isn't already running:
  start();
}

//----------------------------------------------------------------
// the_thread_interface_t::stop_and_go
// 
void
the_thread_interface_t::stop_and_go(std::list<the_transaction_t *> & schedule)
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  
  // remove any further pending transactions:
  while (!transactions_.empty())
  {
    the_transaction_t * t = remove_head(transactions_);
#ifdef DEBUG_THREAD
    cerr << this << " stop_and_go " << t << endl;
#endif
    t->notify(this, the_transaction_t::SKIPPED_E);
  }
  
  // terminate the currently executing transaction:
  terminators().terminate();
  
  // schedule the new transactions:
  transactions_.splice(transactions_.end(), schedule);
  
  // start the thread if it isn't already running:
  start();
}

//----------------------------------------------------------------
// the_thread_interface_t::flush_and_go
// 
void
the_thread_interface_t::flush_and_go(the_transaction_t * transaction)
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  
  // remove any further pending transactions:
  while (!transactions_.empty())
  {
    the_transaction_t * t = remove_head(transactions_);
#ifdef DEBUG_THREAD
    cerr << this << " flush_and_go (1) " << t << endl;
#endif
    t->notify(this, the_transaction_t::SKIPPED_E);
  }
  
  // schedule the next transaction:
  transactions_.push_back(transaction);
  
  // start the thread if it isn't already running:
  start();
}

//----------------------------------------------------------------
// the_thread_interface_t::flush_and_go
// 
void
the_thread_interface_t::
flush_and_go(std::list<the_transaction_t *> & schedule)
{
  the_lock_t<the_mutex_interface_t> locker(mutex_);
  
  // remove any further pending transactions:
  while (!transactions_.empty())
  {
    the_transaction_t * t = remove_head(transactions_);
#ifdef DEBUG_THREAD
    cerr << this << " flush_and_go " << t << endl;
#endif
    t->notify(this, the_transaction_t::SKIPPED_E);
  }
  
  // schedule the next transaction:
  transactions_.splice(transactions_.end(), schedule);
  
  // start the thread if it isn't already running:
  start();
}

//----------------------------------------------------------------
// the_thread_interface_t::work
// 
bool
the_thread_interface_t::work()
{
  the_lock_t<the_mutex_interface_t> lock_this(mutex_, false);
  the_lock_t<the_mutex_interface_t> lock_pool(thread_pool_ ?
					      thread_pool_->mutex_ : nullptr,
					      false);
  bool all_transactions_completed = true;
  
  while (!stopped_)
  {
    // get the next transaction:
    the_transaction_t * t = nullptr;
    {
      lock_pool.arm();
      lock_this.arm();
      
      if (thread_pool_ != nullptr)
      {
	// call back the thread pool:
#ifdef DEBUG_THREAD
	cerr << "thread " << this << " calling the pool" << endl;
#endif
	thread_pool_->handle_thread(thread_pool_cb_data_);
      }
      
      if (transactions_.empty())
      {
	if (sleep_when_idle_ && !stopped_)
	{
	  lock_this.disarm();
	  lock_pool.disarm();
	  take_a_nap(sleep_microsec_);
#ifdef DEBUG_THREAD
	  cerr << this << " sleeping" << endl;
#endif
	  continue;
	}
	else
	{
	  // NOTE: the mutex will remain locked until the function returns:
#ifdef DEBUG_THREAD
	  cerr << "thread " << this << " is finishing up" << endl;
#endif
	  break;
	}
      }
      else
      {
	t = remove_head(transactions_);
#ifdef DEBUG_THREAD
	cerr << "thread " << this << " received " << t << endl;
#endif
	lock_this.disarm();
	lock_pool.disarm();
      }
    }
    
    try
    {
      active_transaction_ = t;
      t->notify(this, the_transaction_t::STARTED_E);
      t->execute(this);
      all_transactions_completed = (all_transactions_completed && true);
      active_transaction_ = nullptr;
      t->notify(this, the_transaction_t::DONE_E);
    }
    catch (std::exception & e)
    {
      active_transaction_ = nullptr;
#ifdef DEBUG_THREAD
      cerr << "FIXME: caught exception: " << e.what() << endl;
#endif
      t->notify(this,
		the_transaction_t::ABORTED_E,
		e.what());
    }
    catch (...)
    {
      active_transaction_ = nullptr;
#ifdef DEBUG_THREAD
      cerr << "FIXME: caught unknonwn exception" << endl;
#endif
      t->notify(this,
		the_transaction_t::ABORTED_E,
		"unknown exception intercepted");
    }
  }
  
  stopped_ = true;
  
  // make sure that all terminators have executed:
#ifndef NDEBUG
  bool ok = the_terminator_t::verify_termination();
  assert(ok);
#endif
  
  all_transactions_completed = (all_transactions_completed &&
				transactions_.empty());
  
  // abort pending transaction:
  while (!transactions_.empty())
  {
    the_transaction_t * t = remove_head(transactions_);
#ifdef DEBUG_THREAD
    cerr << this << " work " << t << endl;
#endif
    t->notify(this, the_transaction_t::SKIPPED_E);
  }
  
#ifdef DEBUG_THREAD
  cerr << "thread " << this << " is finished" << endl;
#endif
  
  if (thread_pool_ != nullptr)
  {
    lock_pool.arm();
    lock_this.arm();
    thread_pool_->handle_thread(thread_pool_cb_data_);
  }
  
  return all_transactions_completed;
}

//----------------------------------------------------------------
// the_thread_interface_t::handle
// 
void
the_thread_interface_t::handle(the_transaction_t * transaction,
			       the_transaction_t::state_t s)
{
  switch (s)
  {
    case the_transaction_t::SKIPPED_E:
    case the_transaction_t::ABORTED_E:
    case the_transaction_t::DONE_E:
      delete transaction;
      break;
      
    default:
      break;
  }
}

//----------------------------------------------------------------
// the_thread_interface_t::blab
// 
void
the_thread_interface_t::blab(const char * message) const
{
  if (thread_pool_ == nullptr)
  {
    cerr << message << endl;
  }
  else
  {
    thread_pool_->blab(message);
  }
}

